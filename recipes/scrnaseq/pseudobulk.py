"""Recipe: pseudobulk aggregation per sample and cell type, and condition DE per cell type with pydeseq2.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §10 (differential expression between conditions); DE steps follow
``RNASEQ_BULK_GUIDELINES.md`` §5-8 through ``recipes/bulk/deseq2.py``.

Procedure
---------
1. ``aggregate_counts``: sum raw counts per (sample, cell type); record the number of cells of each pseudobulk sample.
2. Drop pseudobulk samples with fewer than ``min_cells`` cells (default ``MIN_CELLS``) and report them.
3. ``pseudobulk_de``: per cell type, design-based gene filter and pydeseq2 with the project design; cell types with
   fewer than ``MIN_SAMPLES_PER_GROUP`` samples in any group are skipped and reported.
   Adjusted p-values are given per cell type (``p_adjusted``) and across all cell types (``p_adjusted_global``);
   the project states which one it reports (``STATISTICS_GUIDELINES.md`` §7, multiple testing).

Sample-level covariates (condition, batch, sex, ...) come from ``obs`` and must be constant within a sample.
"""

import logging

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

from recipes.bulk.deseq2 import filter_genes, run_deseq2

logger = logging.getLogger(__name__)

MIN_CELLS = 10  # guide: 10-20 cells per pseudobulk sample
MIN_SAMPLES_PER_GROUP = 2


def aggregate_counts(
    adata: ad.AnnData,
    sample_key: str,
    celltype_key: str,
    covariates: list[str],
    layer: str = "counts",
    min_cells: int = MIN_CELLS,
) -> tuple[dict[str, tuple[pd.DataFrame, pd.DataFrame]], pl.DataFrame]:
    """Sum raw counts per sample and cell type.

    Parameters
    ----------
    adata : ad.AnnData
        Cells after QC and doublet removal, raw counts in ``layer``.
    sample_key : str
        ``obs`` column with the biological sample (the replicate unit).
    celltype_key : str
        ``obs`` column with the final annotation.
    covariates : list[str]
        Sample-level ``obs`` columns carried into the metadata (must be constant within sample).
    layer : str
        Raw counts layer.
    min_cells : int
        Minimum cells per pseudobulk sample.

    Returns
    -------
    pseudobulk : dict[str, tuple[pd.DataFrame, pd.DataFrame]]
        Cell type -> (counts samples x genes, metadata samples x [covariates..., ``n_cells``]).
    report : pl.DataFrame
        One row per (cell type, sample): ``cell_type``, ``sample``, ``n_cells``, ``kept``.

    Raises
    ------
    ValueError
        If a covariate is not constant within a sample.
    """
    obs = adata.obs[[sample_key, celltype_key, *covariates]].astype(str)
    sample_meta = obs.groupby(sample_key, observed=True)[covariates].nunique()
    inconsistent = sample_meta.columns[(sample_meta > 1).any()].to_list()
    if inconsistent:
        raise ValueError(f"covariates not constant within samples: {inconsistent}")
    sample_meta = obs.drop_duplicates(sample_key).set_index(sample_key)[covariates]
    counts = adata.layers[layer]
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)

    pseudobulk, rows = {}, []
    for cell_type, cell_obs in obs.groupby(celltype_key, observed=True):
        samples = cell_obs[sample_key].to_numpy()
        positions = adata.obs_names.get_indexer(cell_obs.index)
        unique_samples, inverse = np.unique(samples, return_inverse=True)
        indicator = sp.csr_matrix(
            (np.ones(samples.size), (inverse, np.arange(samples.size))), shape=(unique_samples.size, samples.size)
        )
        summed = np.asarray((indicator @ counts[positions]).todense())
        n_cells = np.bincount(inverse)
        kept = n_cells >= min_cells
        rows += [
            {"cell_type": cell_type, "sample": s, "n_cells": int(n), "kept": bool(k)}
            for s, n, k in zip(unique_samples, n_cells, kept, strict=True)
        ]
        if not kept.any():
            continue
        frame = pd.DataFrame(summed[kept], index=unique_samples[kept], columns=adata.var_names)
        meta = sample_meta.loc[unique_samples[kept]].copy()
        meta["n_cells"] = n_cells[kept]
        pseudobulk[cell_type] = (frame, meta)
    return pseudobulk, pl.DataFrame(rows)


def pseudobulk_de(
    pseudobulk: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    design: str,
    contrast: tuple[str, str, str],
    lfc_null: float = 0.0,
    min_samples_per_group: int = MIN_SAMPLES_PER_GROUP,
    n_cpus: int = 1,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Condition DE per cell type.

    Returns
    -------
    results : pl.DataFrame
        ``recipes.bulk.deseq2.RESULT_COLUMNS`` plus ``cell_type`` and ``p_adjusted_global`` (BH over all tested genes
        of all cell types).
    skipped : pl.DataFrame
        Cell types not tested and the reason.
    """
    variable = contrast[0]
    tables, skipped = [], []
    for cell_type, (counts, meta) in pseudobulk.items():
        group_sizes = meta[variable].value_counts()
        if set(contrast[1:]) - set(group_sizes.index) or group_sizes.min() < min_samples_per_group:
            skipped.append({"cell_type": cell_type, "reason": f"group sizes {group_sizes.to_dict()}"})
            continue
        genes = filter_genes(counts, meta, variable)
        results, _ = run_deseq2(counts[genes], meta, design, contrast, lfc_null=lfc_null, n_cpus=n_cpus)
        tables.append(results.with_columns(cell_type=pl.lit(cell_type)))
    if not tables:
        return pl.DataFrame(), pl.DataFrame(skipped)
    results = pl.concat(tables)
    p = results["p_value"].to_numpy()
    p_global = np.full(p.size, np.nan)
    finite = np.isfinite(p)
    p_global[finite] = multipletests(p[finite], method="fdr_bh")[1]
    results = results.with_columns(p_adjusted_global=pl.Series(p_global))
    return results, pl.DataFrame(skipped, schema={"cell_type": pl.Utf8, "reason": pl.Utf8})
