"""Recipe: design-based gene filtering, pydeseq2 differential expression with a standard results table, and VST.

Guidelines: ``RNASEQ_BULK_GUIDELINES.md`` §5 (filtering), §6 (VST), §8 (differential expression); used by
``recipes/scrnaseq/pseudobulk.py`` for ``SCRNASEQ_GUIDELINES.md`` §10.

Procedure
---------
1. ``filter_genes``: keep genes with at least ``min_count`` counts in at least as many samples as the smallest group of
   ``group_col``. The rule depends on group sizes only, not on which samples belong to which condition.
2. ``run_deseq2``: ``pydeseq2`` with the given design; Wald test for one contrast, optionally against a minimum effect
   ``lfc_null`` (``alt_hypothesis="greaterAbs"``, the DESeq2 ``lfcThreshold`` test). Returns unshrunken and shrunken
   log fold changes in the column names of the guide.
3. ``vst``: variance-stabilized matrix for exploration, ``use_design=False`` by default (QC).

Limits of ``pydeseq2`` 0.5
--------------------------
No likelihood-ratio test (use DESeq2 in R for factors with several levels or time courses), no gene-level offsets
(import Salmon with ``lengthScaledTPM``, ``recipes/bulk/tximport_lengthscaled.R``).
"""

import logging

import numpy as np
import pandas as pd
import polars as pl
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

logger = logging.getLogger(__name__)

MIN_COUNT = 10
RESULT_COLUMNS = [
    "gene_id",
    "base_mean",
    "log2_fold_change",
    "lfc_se",
    "log2_fold_change_shrunk",
    "stat",
    "p_value",
    "p_adjusted",
    "contrast",
]


def filter_genes(counts: pd.DataFrame, metadata: pd.DataFrame, group_col: str, min_count: int = MIN_COUNT) -> pd.Index:
    """Return genes passing the design-based filter.

    Parameters
    ----------
    counts : pd.DataFrame
        Samples x genes, raw counts.
    metadata : pd.DataFrame
        Samples x variables, same index as ``counts``.
    group_col : str
        Column defining the groups of the design (e.g. condition).
    min_count : int
        Minimum count per sample.

    Returns
    -------
    pd.Index
        Genes kept.
    """
    if not counts.index.equals(metadata.index):
        raise ValueError("counts and metadata must have the same sample index, in the same order")
    n_required = int(metadata[group_col].value_counts().min())
    keep = (counts >= min_count).sum(axis=0) >= n_required
    logger.info(
        f"gene filter: {int(keep.sum())} of {keep.size} genes (>= {min_count} counts in >= {n_required} samples)"
    )
    return counts.columns[keep.to_numpy()]


def run_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    contrast: tuple[str, str, str],
    lfc_null: float = 0.0,
    alpha: float = 0.05,
    n_cpus: int = 1,
    shrink: bool = True,
) -> tuple[pl.DataFrame, DeseqDataSet]:
    """Fit pydeseq2 and return the standard results table for one Wald contrast.

    Parameters
    ----------
    counts : pd.DataFrame
        Samples x genes, raw (or ``lengthScaledTPM``) integer counts, already filtered.
    metadata : pd.DataFrame
        Samples x covariates; categorical variables as strings.
    design : str
        Formula, e.g. ``"~ batch + sex + condition"``; the variable of interest last.
    contrast : tuple[str, str, str]
        ``(variable, tested_level, reference_level)``.
    lfc_null : float
        Minimum absolute log2 fold change tested (0 tests for any change).
    alpha : float
        FDR level used by independent filtering.
    n_cpus : int
        CPUs for pydeseq2.
    shrink : bool
        Whether to add apeGLM-type shrunken log fold changes.

    Returns
    -------
    results : pl.DataFrame
        Columns ``RESULT_COLUMNS`` (``gene_symbol`` is added by the caller from the annotation).
    dds : DeseqDataSet
        Fitted object (for VST or further contrasts).
    """
    variable, tested, reference = contrast
    integer_counts = counts.round().astype(int)
    metadata = metadata.copy()
    levels = [reference] + sorted(set(metadata[variable].astype(str)) - {reference})
    metadata[variable] = pd.Categorical(metadata[variable].astype(str), categories=levels)  # reference level first
    dds = DeseqDataSet(counts=integer_counts, metadata=metadata, design=design, n_cpus=n_cpus, quiet=True)
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=[variable, tested, reference],
        alpha=alpha,
        lfc_null=lfc_null,
        alt_hypothesis="greaterAbs" if lfc_null > 0 else None,
        n_cpus=n_cpus,
        quiet=True,
    )
    stats.summary()
    raw = stats.results_df.copy()
    shrunk = pd.Series(np.nan, index=raw.index)
    if shrink:
        coeff = _coefficient_name(dds, variable, tested, reference)
        stats.lfc_shrink(coeff=coeff)
        shrunk = stats.results_df["log2FoldChange"]
    results = pl.DataFrame(
        {
            "gene_id": raw.index.to_list(),
            "base_mean": raw["baseMean"].to_numpy(),
            "log2_fold_change": raw["log2FoldChange"].to_numpy(),
            "lfc_se": raw["lfcSE"].to_numpy(),
            "log2_fold_change_shrunk": shrunk.reindex(raw.index).to_numpy(),
            "stat": raw["stat"].to_numpy(),
            "p_value": raw["pvalue"].to_numpy(),
            "p_adjusted": raw["padj"].to_numpy(),
            "contrast": f"{variable}_{tested}_vs_{reference}",
        }
    )
    return results.select(RESULT_COLUMNS), dds


def _coefficient_name(dds: DeseqDataSet, variable: str, tested: str, reference: str) -> str:
    """Name of the design-matrix column for ``tested`` vs ``reference`` (formulaic naming)."""
    candidates = [c for c in dds.obsm["design_matrix"].columns if c.startswith(f"{variable}[") and tested in c]
    if len(candidates) != 1:
        raise ValueError(
            f"cannot find a unique coefficient for {variable} {tested} vs {reference}: {candidates}. Shrinkage needs"
            f" {reference!r} as the reference level of {variable!r} in the design."
        )
    return candidates[0]


def vst(dds: DeseqDataSet, use_design: bool = False) -> pd.DataFrame:
    """Variance-stabilized expression (samples x genes) from a pydeseq2 object."""
    dds.vst(use_design=use_design)
    return pd.DataFrame(dds.layers["vst_counts"], index=dds.obs_names, columns=dds.var_names)
