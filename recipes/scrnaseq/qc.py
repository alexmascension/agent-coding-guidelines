"""Recipe: gene-group flags, nuclear fraction and adaptive per-capture cell filters.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §4.1 (cell calling, nuclear fraction) and §4.2 (cell-level
filters).

Procedure
---------
1. ``flag_gene_groups``: ``var["mt"]`` from the chromosome in the annotation (never from a symbol
   prefix);
   ``var["ribo"]`` and ``var["hb"]`` from explicit species-specific ID lists stored in the project.
   Targeted panels
   without mitochondrial genes pass ``require_mito=False``.
2. ``add_nuclear_fraction``: ``obs["nuclear_fraction"] = unspliced / (spliced + unspliced)`` from
   the STARsolo
   ``Velocyto`` layers.
3. ``qc_metrics`` on the raw ``counts`` layer (``SCRNASEQ_GUIDELINES.md`` §1.1 states the layer of
   each step).
4. ``mad_outlier_flags``: per capture, flags cells more than ``n`` scaled MADs from the capture
   median. The MAD is
   scaled by 1.4826 so that it estimates the SD under normality: 5 scaled MADs = 7.41 unscaled
   MADs. Directions: both
   tails for log total counts and log genes detected; upper tail for mitochondrial fraction; lower
   tail for nuclear
   fraction.

Policies that the caller sees in the report
-------------------------------------------
* **Non-finite values** (e.g. nuclear fraction of a barcode with no spliced or unspliced reads)
  cannot be evaluated:
  they are flagged ``qc_invalid_<metric>``, never counted as passing. ``qc_fail = qc_outlier |
  qc_invalid``.
* **MAD = 0** (typical for mitochondrial fraction in nuclei, where most barcodes have none): the
  MAD rule cannot be
  applied. The absolute threshold ``max_mito_fraction`` is used if given (project value justified
  by tissue and
  assay); otherwise the metric is reported as ``not_evaluated_mad_zero`` and a warning is logged.
* ``max_mito_fraction``, when given, is applied in every capture in addition to the MAD rule.

Defaults
--------
``N_MADS_COUNTS = 5`` and ``N_MADS_MITO = 3`` are within the 3-5 range of the guide; projects state
their own values.
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
REPORT_SCHEMA = {
    "capture": pl.Utf8,
    "metric": pl.Utf8,
    "median": pl.Float64,
    "scaled_mad": pl.Float64,
    "lower": pl.Float64,
    "upper": pl.Float64,
    "status": pl.Utf8,
    "n_flagged": pl.Int64,
    "n_invalid": pl.Int64,
    "n_cells": pl.Int64,
}


def flag_gene_groups(
    adata: ad.AnnData,
    annotation: pl.DataFrame,
    ribosomal_ids: list[str] | None = None,
    haemoglobin_ids: list[str] | None = None,
    require_mito: bool = True,
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
    require_mito : bool
        Raise if no mitochondrial gene is present (whole-transcriptome data); False for targeted
        panels.

    Returns
    -------
    dict[str, int]
        Number of genes flagged per group (to record in the QC report).

    Raises
    ------
    ValueError
        If ``require_mito`` and no mitochondrial gene is found (wrong annotation or chromosome
        naming).
    """
    mito_ids = set(
        annotation.filter(pl.col("chromosome").is_in(MITOCHONDRIAL_CHROMOSOMES))[
            "gene_id"
        ].to_list()
    )
    adata.var["mt"] = adata.var_names.isin(mito_ids)
    adata.var["ribo"] = adata.var_names.isin(set(ribosomal_ids or []))
    adata.var["hb"] = adata.var_names.isin(set(haemoglobin_ids or []))
    counts = {group: int(adata.var[group].sum()) for group in ("mt", "ribo", "hb")}
    if counts["mt"] == 0 and require_mito:
        raise ValueError(
            f"no mitochondrial genes found on chromosomes {MITOCHONDRIAL_CHROMOSOMES}"
        )
    logger.info(f"gene groups flagged: {counts}")
    return counts


def add_nuclear_fraction(
    adata: ad.AnnData, spliced: str = "spliced", unspliced: str = "unspliced"
) -> None:
    """Add ``obs["nuclear_fraction"]`` = unspliced / (spliced + unspliced); NaN if both are 0."""
    s = np.asarray(adata.layers[spliced].sum(axis=1)).ravel().astype(float)
    u = np.asarray(adata.layers[unspliced].sum(axis=1)).ravel().astype(float)
    total = s + u
    with np.errstate(invalid="ignore", divide="ignore"):
        adata.obs["nuclear_fraction"] = np.where(total > 0, u / total, np.nan)


def qc_metrics(adata: ad.AnnData, layer: str = "counts") -> None:
    """Add ``obs`` QC columns (log1p totals and genes; ``mito_fraction`` if mito genes exist)."""
    x = adata.layers[layer]
    x = x.tocsr() if sp.issparse(x) else np.asarray(x)
    total = np.asarray(x.sum(axis=1)).ravel().astype(float)
    n_genes = np.asarray((x > 0).sum(axis=1)).ravel().astype(float)
    adata.obs["log1p_total_counts"] = np.log1p(total)
    adata.obs["log1p_n_genes"] = np.log1p(n_genes)
    if "mt" in adata.var and adata.var["mt"].any():
        mito = (
            np.asarray(x[:, np.flatnonzero(adata.var["mt"].to_numpy())].sum(axis=1))
            .ravel()
            .astype(float)
        )
        with np.errstate(invalid="ignore", divide="ignore"):
            adata.obs["mito_fraction"] = np.where(total > 0, mito / total, np.nan)


def _median_scaled_mad(values: np.ndarray) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan"), float("nan")
    median = float(np.median(finite))
    return median, MAD_SCALE * float(np.median(np.abs(finite - median)))


def mad_outlier_flags(
    adata: ad.AnnData,
    capture_key: str,
    n_mads_counts: float = N_MADS_COUNTS,
    n_mads_mito: float = N_MADS_MITO,
    n_mads_nuclear: float | None = N_MADS_NUCLEAR,
    max_mito_fraction: float | None = None,
) -> pl.DataFrame:
    """Flag per-capture outliers and non-evaluable values (in place).

    Adds ``obs["qc_outlier_<metric>"]``, ``obs["qc_invalid_<metric>"]``, ``obs["qc_outlier"]``,
    ``obs["qc_invalid"]``
    and ``obs["qc_fail"]``. Requires ``qc_metrics`` (and ``add_nuclear_fraction`` unless
    ``n_mads_nuclear`` is None).
    The mitochondrial rule is skipped when ``mito_fraction`` is absent (panels without
    mitochondrial genes).

    Parameters
    ----------
    adata : ad.AnnData
        Cells of one or more captures.
    capture_key : str
        ``obs`` column identifying the capture (library); thresholds are computed within each
        capture.
    n_mads_counts, n_mads_mito, n_mads_nuclear : float
        Width of the accepted interval in scaled MADs; ``n_mads_nuclear=None`` skips the
        nuclear-fraction filter.
    max_mito_fraction : float | None
        Absolute upper limit for the mitochondrial fraction, applied in every capture (project
        value).

    Returns
    -------
    pl.DataFrame
        One row per capture x metric (schema ``REPORT_SCHEMA``); ``status`` is ``evaluated``,
        ``absolute_threshold_only`` or ``not_evaluated_mad_zero``.
    """
    rules = [
        ("log1p_total_counts", n_mads_counts, True, True),
        ("log1p_n_genes", n_mads_counts, True, True),
    ]
    if "mito_fraction" in adata.obs:
        rules.append(("mito_fraction", n_mads_mito, False, True))
    if n_mads_nuclear is not None:
        rules.append(("nuclear_fraction", n_mads_nuclear, True, False))
    outlier = {metric: np.zeros(adata.n_obs, dtype=bool) for metric, *_ in rules}
    invalid = {metric: np.zeros(adata.n_obs, dtype=bool) for metric, *_ in rules}
    captures = adata.obs[capture_key].astype(str).to_numpy()
    rows = []
    for capture in np.unique(captures):
        in_capture = captures == capture
        for metric, n_mads, use_lower, use_upper in rules:
            values = adata.obs[metric].to_numpy(dtype=float)[in_capture]
            finite = np.isfinite(values)
            median, mad = _median_scaled_mad(values)
            lower = median - n_mads * mad if use_lower else None
            upper = median + n_mads * mad if use_upper else None
            status = "evaluated"
            if not mad > 0:
                lower = upper = None
                status = "not_evaluated_mad_zero"
            flagged = np.zeros(values.size, dtype=bool)
            if lower is not None:
                flagged |= finite & (values < lower)
            if upper is not None:
                flagged |= finite & (values > upper)
            if metric == "mito_fraction" and max_mito_fraction is not None:
                flagged |= finite & (values > max_mito_fraction)
                if status == "not_evaluated_mad_zero":
                    status = "absolute_threshold_only"
                upper = max_mito_fraction if upper is None else min(upper, max_mito_fraction)
            if status == "not_evaluated_mad_zero":
                logger.warning(f"capture {capture}: MAD of {metric} is 0; metric not evaluated")
            outlier[metric][in_capture] = flagged
            invalid[metric][in_capture] = ~finite
            rows.append(
                {
                    "capture": capture,
                    "metric": metric,
                    "median": median,
                    "scaled_mad": mad,
                    "lower": lower,
                    "upper": upper,
                    "status": status,
                    "n_flagged": int(flagged.sum()),
                    "n_invalid": int((~finite).sum()),
                    "n_cells": int(in_capture.sum()),
                }
            )
    any_outlier = np.zeros(adata.n_obs, dtype=bool)
    any_invalid = np.zeros(adata.n_obs, dtype=bool)
    for metric in outlier:
        adata.obs[f"qc_outlier_{metric}"] = outlier[metric]
        adata.obs[f"qc_invalid_{metric}"] = invalid[metric]
        any_outlier |= outlier[metric]
        any_invalid |= invalid[metric]
    adata.obs["qc_outlier"] = any_outlier
    adata.obs["qc_invalid"] = any_invalid
    adata.obs["qc_fail"] = any_outlier | any_invalid
    logger.info(
        f"QC: {int(any_outlier.sum())} outliers, {int(any_invalid.sum())} invalid, "
        f"of {adata.n_obs} cells"
    )
    return pl.DataFrame(rows, schema=REPORT_SCHEMA)
