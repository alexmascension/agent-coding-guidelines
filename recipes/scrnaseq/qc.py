"""Recipe: gene-group flags, nuclear fraction and adaptive per-capture cell filters.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §4.1 (cell calling, nuclear fraction) and §4.2 (cell-level filters).

Procedure
---------
1. ``flag_gene_groups``: ``var["mt"]`` from the chromosome in the annotation (never from a symbol prefix);
   ``var["ribo"]`` and ``var["hb"]`` from explicit species-specific ID lists stored in the project.
2. ``add_nuclear_fraction``: ``obs["nuclear_fraction"] = unspliced / (spliced + unspliced)`` from the STARsolo
   ``Velocyto`` layers.
3. ``mad_outlier_flags``: per capture, flags cells more than ``n_mads`` MADs (median absolute deviations, scaled to be
   consistent with the SD under normality) from the capture median. Directions: both tails for log total counts and
   log genes detected; upper tail for mitochondrial fraction; lower tail for nuclear fraction.
   Returns the per-capture thresholds and the number of cells flagged by each metric.

Defaults
--------
``N_MADS_COUNTS = 5`` and ``N_MADS_MITO = 3`` are within the 3-5 range of the guide; projects state their own values.
The mitochondrial filter is applied only when the capture median is above zero (nuclei data have ~0 mitochondrial
reads, where the MAD is 0 and any read would be an outlier).
"""

import logging

import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp

from recipes.common.gene_ids import MITOCHONDRIAL_CHROMOSOMES

logger = logging.getLogger(__name__)

N_MADS_COUNTS = 5
N_MADS_MITO = 3
N_MADS_NUCLEAR = 3
MAD_SCALE = 1.4826  # makes the MAD a consistent estimator of the SD for normal data


def flag_gene_groups(
    adata: ad.AnnData,
    annotation: "pl.DataFrame",
    ribosomal_ids: list[str] | None = None,
    haemoglobin_ids: list[str] | None = None,
) -> dict[str, int]:
    """Add boolean ``var`` columns ``mt``, ``ribo`` and ``hb`` (in place).

    Parameters
    ----------
    adata : ad.AnnData
        ``var_names`` are stable Ensembl IDs (``recipes.common.gene_ids.set_var_gene_ids``).
    annotation : pl.DataFrame
        Annotation table with ``gene_id`` and ``chromosome``.
    ribosomal_ids, haemoglobin_ids : list[str] | None
        Explicit ID lists stored in the project; None leaves the column all False.

    Returns
    -------
    dict[str, int]
        Number of genes flagged per group (to record in the QC report).

    Raises
    ------
    ValueError
        If no mitochondrial gene is found (wrong annotation or chromosome naming).
    """
    mito_ids = set(annotation.filter(pl.col("chromosome").is_in(MITOCHONDRIAL_CHROMOSOMES))["gene_id"].to_list())
    adata.var["mt"] = adata.var_names.isin(mito_ids)
    adata.var["ribo"] = adata.var_names.isin(set(ribosomal_ids or []))
    adata.var["hb"] = adata.var_names.isin(set(haemoglobin_ids or []))
    counts = {group: int(adata.var[group].sum()) for group in ("mt", "ribo", "hb")}
    if counts["mt"] == 0:
        raise ValueError(f"no mitochondrial genes found on chromosomes {MITOCHONDRIAL_CHROMOSOMES}")
    logger.info(f"gene groups flagged: {counts}")
    return counts


def add_nuclear_fraction(adata: ad.AnnData, spliced: str = "spliced", unspliced: str = "unspliced") -> None:
    """Add ``obs["nuclear_fraction"]`` = unspliced / (spliced + unspliced) per barcode (NaN if both are 0)."""
    s = np.asarray(adata.layers[spliced].sum(axis=1)).ravel().astype(float)
    u = np.asarray(adata.layers[unspliced].sum(axis=1)).ravel().astype(float)
    total = s + u
    with np.errstate(invalid="ignore", divide="ignore"):
        adata.obs["nuclear_fraction"] = np.where(total > 0, u / total, np.nan)


def qc_metrics(adata: ad.AnnData, layer: str = "counts") -> None:
    """Add ``obs`` columns ``log1p_total_counts``, ``log1p_n_genes``, ``mito_fraction`` from raw counts (in place)."""
    x = adata.layers[layer]
    x = x.tocsr() if sp.issparse(x) else np.asarray(x)
    total = np.asarray(x.sum(axis=1)).ravel().astype(float)
    n_genes = np.asarray((x > 0).sum(axis=1)).ravel().astype(float)
    mito = np.asarray(x[:, np.flatnonzero(adata.var["mt"].to_numpy())].sum(axis=1)).ravel().astype(float)
    adata.obs["log1p_total_counts"] = np.log1p(total)
    adata.obs["log1p_n_genes"] = np.log1p(n_genes)
    with np.errstate(invalid="ignore", divide="ignore"):
        adata.obs["mito_fraction"] = np.where(total > 0, mito / total, np.nan)


def _bounds(values: np.ndarray, n_mads: float) -> tuple[float, float, float]:
    finite = values[np.isfinite(values)]
    median = float(np.median(finite))
    mad = MAD_SCALE * float(np.median(np.abs(finite - median)))
    return median, median - n_mads * mad, median + n_mads * mad


def mad_outlier_flags(
    adata: ad.AnnData,
    capture_key: str,
    n_mads_counts: float = N_MADS_COUNTS,
    n_mads_mito: float = N_MADS_MITO,
    n_mads_nuclear: float | None = N_MADS_NUCLEAR,
) -> pl.DataFrame:
    """Flag per-capture outliers (in place: ``obs["qc_outlier_<metric>"]`` and ``obs["qc_outlier"]``).

    Requires ``qc_metrics`` (and ``add_nuclear_fraction`` unless ``n_mads_nuclear`` is None).

    Parameters
    ----------
    adata : ad.AnnData
        Cells of one or more captures.
    capture_key : str
        ``obs`` column identifying the capture (library); thresholds are computed within each capture.
    n_mads_counts, n_mads_mito, n_mads_nuclear : float
        Width of the accepted interval in scaled MADs; ``n_mads_nuclear=None`` skips the nuclear-fraction filter.

    Returns
    -------
    pl.DataFrame
        One row per capture x metric: ``capture``, ``metric``, ``median``, ``lower``, ``upper`` (null when that tail is
        not used), ``n_flagged``, ``n_cells``.
    """
    rules = [
        ("log1p_total_counts", n_mads_counts, True, True),
        ("log1p_n_genes", n_mads_counts, True, True),
        ("mito_fraction", n_mads_mito, False, True),
    ]
    if n_mads_nuclear is not None:
        rules.append(("nuclear_fraction", n_mads_nuclear, True, False))
    flags = {metric: np.zeros(adata.n_obs, dtype=bool) for metric, *_ in rules}
    captures = adata.obs[capture_key].astype(str).to_numpy()
    rows = []
    for capture in np.unique(captures):
        in_capture = captures == capture
        for metric, n_mads, use_lower, use_upper in rules:
            values = adata.obs[metric].to_numpy(dtype=float)[in_capture]
            median, lower, upper = _bounds(values, n_mads)
            if metric == "mito_fraction" and median == 0:
                use_upper = False  # nuclei: MAD is 0, every mitochondrial read would be flagged
            flagged = np.zeros(values.size, dtype=bool)
            if use_lower:
                flagged |= values < lower
            if use_upper:
                flagged |= values > upper
            flags[metric][in_capture] = flagged
            rows.append(
                {
                    "capture": capture,
                    "metric": metric,
                    "median": median,
                    "lower": lower if use_lower else None,
                    "upper": upper if use_upper else None,
                    "n_flagged": int(flagged.sum()),
                    "n_cells": int(in_capture.sum()),
                }
            )
    any_flag = np.zeros(adata.n_obs, dtype=bool)
    for metric, values in flags.items():
        adata.obs[f"qc_outlier_{metric}"] = values
        any_flag |= values
    adata.obs["qc_outlier"] = any_flag
    logger.info(f"QC outliers: {int(any_flag.sum())} of {adata.n_obs} cells")
    return pl.DataFrame(
        rows,
        schema={
            "capture": pl.Utf8,
            "metric": pl.Utf8,
            "median": pl.Float64,
            "lower": pl.Float64,
            "upper": pl.Float64,
            "n_flagged": pl.Int64,
            "n_cells": pl.Int64,
        },
    )
