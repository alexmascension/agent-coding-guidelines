"""Recipe: pseudobulk aggregation per unit and cell type, and condition DE per cell type.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §10 (differential expression between conditions); DE steps
follow
``RNASEQ_BULK_GUIDELINES.md`` §5-8 through ``recipes/bulk/deseq2.py``.

Procedure
---------
1. ``aggregate_counts``: sum counts per (aggregation unit, cell type) and record the number of
   cells. The aggregation
   unit is the experimental unit of the comparison: a sample, or a donor only when each donor
   contributes to one
   level of every compared factor. When a donor has samples in several conditions or time points,
   aggregate per
   sample and model the donor (``SCRNASEQ_GUIDELINES.md`` §2).
2. Units with fewer than ``min_cells`` cells (default ``MIN_CELLS``) are not used. The report lists
   every
   (unit, cell type) combination with status ``kept``, ``below_min_cells`` or ``absent`` (no
   cells). A cell type
   absent from a unit is not turned into a pseudobulk sample of zeros.
3. ``pseudobulk_de``: per cell type, design-based gene filter and pydeseq2 with the project design;
   cell types where
   a compared group has fewer than ``MIN_SAMPLES_PER_GROUP`` units are reported as
   ``not_testable``.
   Adjusted p-values are given per cell type (``p_adjusted``, pydeseq2 with independent filtering
   within the cell
   type) and across all cell types (``p_adjusted_global``, BH over the raw p-values of all tested
   genes). The two
   families differ, also because independent filtering is applied per cell type; the project states
   which one it
   reports (``STATISTICS_GUIDELINES.md`` §7, multiple testing).

Covariates keep their dtype: numeric covariates (age, RIN, postmortem interval) stay numeric and
enter the design as
slopes; categorical ones stay categorical. Missing values MUST be resolved before calling (an error
is raised).
"""

import logging

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

from recipes.bulk.deseq2 import RESULT_COLUMNS, filter_genes, run_deseq2

logger = logging.getLogger(__name__)

MIN_CELLS = 10  # guide: 10-20 cells per pseudobulk sample
MIN_SAMPLES_PER_GROUP = (
    2  # computational minimum for a dispersion estimate, not a sample-size justification
)
REPORT_SCHEMA = {"cell_type": pl.Utf8, "unit": pl.Utf8, "n_cells": pl.Int64, "status": pl.Utf8}


def unit_metadata(obs: pd.DataFrame, unit_key: str, covariates: list[str]) -> pd.DataFrame:
    """One row per aggregation unit with its covariates, keeping dtypes.

    Raises
    ------
    ValueError
        If a covariate has missing values or is not constant within a unit.
    """
    missing = [c for c in covariates if obs[c].isna().any()]
    if missing:
        raise ValueError(
            f"covariates with missing values: {missing}; resolve them before aggregation"
        )
    units = obs[unit_key].astype(str)
    varying = [
        c
        for c in covariates
        if (obs.groupby(units, observed=True)[c].nunique(dropna=False) > 1).any()
    ]
    if varying:
        raise ValueError(f"covariates not constant within {unit_key}: {varying}")
    meta = obs.assign(_unit=units).drop_duplicates("_unit").set_index("_unit")[covariates]
    meta.index.name = unit_key
    for column in covariates:
        if isinstance(meta[column].dtype, pd.CategoricalDtype):
            meta[column] = meta[column].cat.remove_unused_categories()
    return meta


def aggregate_counts(
    adata: ad.AnnData,
    unit_key: str,
    celltype_key: str,
    covariates: list[str],
    layer: str = "counts",
    min_cells: int = MIN_CELLS,
) -> tuple[dict[str, tuple[pd.DataFrame, pd.DataFrame]], pl.DataFrame]:
    """Sum counts per aggregation unit and cell type.

    Parameters
    ----------
    adata : ad.AnnData
        Cells after QC and doublet removal; counts in ``layer``.
    unit_key : str
        ``obs`` column with the experimental unit (sample, or donor when appropriate; see module
        docstring).
    celltype_key : str
        ``obs`` column with the final annotation.
    covariates : list[str]
        Unit-level ``obs`` columns carried into the metadata; dtypes are preserved.
    layer : str
        Counts layer (``SCRNASEQ_GUIDELINES.md`` §1.1 states which one each step uses).
    min_cells : int
        Minimum cells per pseudobulk sample.

    Returns
    -------
    pseudobulk : dict[str, tuple[pd.DataFrame, pd.DataFrame]]
        Cell type -> (counts units x genes, metadata units x [covariates..., ``n_cells``]), kept
        units only.
    report : pl.DataFrame
        Every (cell type, unit) combination: ``cell_type``, ``unit``, ``n_cells``, ``status``.
    """
    meta_all = unit_metadata(adata.obs, unit_key, covariates)
    units = adata.obs[unit_key].astype(str).to_numpy()
    cell_types = adata.obs[celltype_key].astype(str).to_numpy()
    counts = adata.layers[layer]
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)
    all_units = meta_all.index.to_numpy()

    pseudobulk, rows = {}, []
    for cell_type in np.unique(cell_types):
        in_type = np.flatnonzero(cell_types == cell_type)
        unit_index = pd.Index(all_units).get_indexer(units[in_type])
        indicator = sp.csr_matrix(
            (np.ones(in_type.size), (unit_index, np.arange(in_type.size))),
            shape=(all_units.size, in_type.size),
        )
        summed = (indicator @ counts[in_type]).toarray()
        n_cells = np.bincount(unit_index, minlength=all_units.size)
        kept = n_cells >= min_cells
        for unit, n, k in zip(all_units, n_cells, kept, strict=True):
            status = "kept" if k else ("absent" if n == 0 else "below_min_cells")
            rows.append(
                {"cell_type": cell_type, "unit": unit, "n_cells": int(n), "status": status}
            )
        if not kept.any():
            continue
        frame = pd.DataFrame(summed[kept], index=all_units[kept], columns=adata.var_names)
        meta = meta_all.loc[all_units[kept]].copy()
        meta["n_cells"] = n_cells[kept]
        pseudobulk[cell_type] = (frame, meta)
    return pseudobulk, pl.DataFrame(rows, schema=REPORT_SCHEMA)


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
        ``recipes.bulk.deseq2.RESULT_COLUMNS`` plus ``cell_type`` and ``p_adjusted_global``; empty
        with the same
        schema when nothing was testable.
    not_testable : pl.DataFrame
        ``cell_type``, ``reason``.
    """
    variable, tested, reference = contrast
    tables, skipped = [], []
    for cell_type, (counts, meta) in pseudobulk.items():
        group_sizes = meta[variable].astype(str).value_counts()
        sizes = {level: int(group_sizes.get(level, 0)) for level in (tested, reference)}
        if min(sizes.values()) < min_samples_per_group:
            skipped.append({"cell_type": cell_type, "reason": f"units per group {sizes}"})
            continue
        genes = filter_genes(counts, meta, variable)
        results, _ = run_deseq2(
            counts[genes], meta, design, contrast, lfc_null=lfc_null, n_cpus=n_cpus
        )
        tables.append(results.with_columns(cell_type=pl.lit(cell_type)))
    skipped_table = pl.DataFrame(skipped, schema={"cell_type": pl.Utf8, "reason": pl.Utf8})
    if not tables:
        empty = pl.DataFrame(
            schema={
                **{c: pl.Float64 for c in RESULT_COLUMNS},
                "gene_id": pl.Utf8,
                "contrast": pl.Utf8,
                "cell_type": pl.Utf8,
                "p_adjusted_global": pl.Float64,
            }
        )
        return empty.select([*RESULT_COLUMNS, "cell_type", "p_adjusted_global"]), skipped_table
    results = pl.concat(tables)
    p = results["p_value"].to_numpy()
    p_global = np.full(p.size, np.nan)
    finite = np.isfinite(p)
    p_global[finite] = multipletests(p[finite], method="fdr_bh")[1]
    return results.with_columns(p_adjusted_global=pl.Series(p_global)), skipped_table
