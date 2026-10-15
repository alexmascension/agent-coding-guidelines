"""Recipe: design-based gene filter, pydeseq2 DE with a standard results table, and VST.

Guidelines: ``RNASEQ_BULK_GUIDELINES.md`` §5 (filtering), §6 (VST), §8 (differential expression);
used by
``recipes/scrnaseq/pseudobulk.py`` for ``SCRNASEQ_GUIDELINES.md`` §10.

Procedure
---------
1. ``filter_genes``: keep genes with at least ``min_count`` counts in at least as many samples as
   the smallest
   **observed** group of ``group_col``. The rule depends on group sizes only.
2. ``run_deseq2``: ``pydeseq2`` with the given design; Wald test for one two-level contrast of a
   main effect,
   optionally against a minimum effect ``lfc_null`` (``alt_hypothesis="greaterAbs"``, the DESeq2
   ``lfcThreshold``
   test). Returns unshrunken and shrunken log fold changes in the column names of the guide.
   Meta-analyses use the
   unshrunken estimates and their standard errors (``META_ANALYSIS_GUIDELINES.md``).
3. ``vst``: variance-stabilized matrix for exploration, ``use_design=False`` by default (QC).

Input contract
--------------
Counts are non-negative and finite. Non-integer values are accepted only with
``allow_non_integer=True``, for
``lengthScaledTPM`` counts from tximport (rounded, as DESeq2 users do); normalized or transformed
values MUST NOT be
passed: rounding does not make them counts.

Limits of ``pydeseq2`` 0.5
--------------------------
No likelihood-ratio test (use DESeq2 in R for factors with several levels or time courses); no
gene-level offsets
(import Salmon with ``lengthScaledTPM``); shrinkage only for a main-effect coefficient of the
treatment-coded design,
so designs with interactions involving the contrasted variable are rejected here.
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


def filter_genes(
    counts: pd.DataFrame, metadata: pd.DataFrame, group_col: str, min_count: int = MIN_COUNT
) -> pd.Index:
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

    Raises
    ------
    ValueError
        If indices differ, metadata is empty, or fewer than two groups are observed.
    """
    if metadata.empty or not counts.index.equals(metadata.index):
        raise ValueError(
            "counts and metadata must be non-empty with the same sample index, in the same order"
        )
    group_sizes = (
        metadata[group_col].astype(str).value_counts()
    )  # observed levels only (unused categories dropped)
    if group_sizes.size < 2:
        raise ValueError(
            f"fewer than two observed groups in {group_col!r}: {group_sizes.to_dict()}"
        )
    n_required = int(group_sizes.min())
    keep = (counts >= min_count).sum(axis=0) >= n_required
    logger.info(
        f"gene filter: {int(keep.sum())}/{keep.size} genes"
        f" (>= {min_count} counts in >= {n_required} samples)"
    )
    return counts.columns[keep.to_numpy()]


def validated_counts(counts: pd.DataFrame, allow_non_integer: bool = False) -> pd.DataFrame:
    """Return integer counts after checking that the input is a count matrix.

    Raises
    ------
    ValueError
        If values are negative or non-finite, or non-integer while ``allow_non_integer`` is False.
    """
    values = counts.to_numpy(dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("counts contain negative or non-finite values")
    if not allow_non_integer and not np.allclose(values, np.round(values)):
        raise ValueError(
            "non-integer counts; pass allow_non_integer=True only for lengthScaledTPM counts"
        )
    return counts.round().astype(np.int64)


def coefficient_name(design_columns: list[str], variable: str, tested: str) -> str:
    """Exact formulaic name of the treatment-coded coefficient ``variable[T.tested]``.

    Raises
    ------
    ValueError
        If the coefficient is absent or the variable takes part in an interaction.
    """
    interactions = [
        c
        for c in design_columns
        if ":" in c and any(part.startswith(f"{variable}[") for part in c.split(":"))
    ]
    if interactions:
        raise ValueError(
            f"{variable!r} is in interactions {interactions}; specify the contrast in R (DESeq2)"
        )
    name = f"{variable}[T.{tested}]"
    if name not in design_columns:
        raise ValueError(f"coefficient {name!r} not in the design matrix columns {design_columns}")
    return name


def run_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    contrast: tuple[str, str, str],
    lfc_null: float = 0.0,
    alpha: float = 0.05,
    n_cpus: int = 1,
    shrink: bool = True,
    allow_non_integer: bool = False,
) -> tuple[pl.DataFrame, DeseqDataSet]:
    """Fit pydeseq2 and return the standard results table for one Wald contrast.

    Parameters
    ----------
    counts : pd.DataFrame
        Samples x genes counts, already filtered (see module docstring for the input contract).
    metadata : pd.DataFrame
        Samples x covariates. Numeric columns enter the design as slopes; string or categorical
        columns as factors.
    design : str
        Formula, e.g. ``"~ batch + sex + age + condition"``.
    contrast : tuple[str, str, str]
        ``(variable, tested_level, reference_level)``; the reference becomes the first level of
        ``variable``.
    lfc_null : float
        Minimum absolute log2 fold change tested (0 tests for any change).
    alpha : float
        FDR level used by independent filtering.
    n_cpus : int
        CPUs for pydeseq2.
    shrink : bool
        Whether to add apeGLM-type shrunken log fold changes.
    allow_non_integer : bool
        Accept non-integer counts (``lengthScaledTPM``) and round them.

    Returns
    -------
    results : pl.DataFrame
        Columns ``RESULT_COLUMNS`` (``gene_symbol`` is added by the caller from the annotation).
    dds : DeseqDataSet
        Fitted object (for VST or further contrasts).

    Raises
    ------
    ValueError
        If the counts violate the input contract, a covariate has missing values, or a level is
        absent.
    """
    variable, tested, reference = contrast
    integer_counts = validated_counts(counts, allow_non_integer)
    metadata = metadata.copy()
    if metadata.isna().any().any():
        raise ValueError(
            f"metadata has missing values in {metadata.columns[metadata.isna().any()].to_list()}"
        )
    present = set(metadata[variable].astype(str))
    if not {tested, reference} <= present:
        raise ValueError(
            f"levels {tested!r} and {reference!r} must be present in {variable!r}; found {present}"
        )
    levels = [reference] + sorted(present - {reference})
    metadata[variable] = pd.Categorical(
        metadata[variable].astype(str), categories=levels
    )  # reference level first
    dds = DeseqDataSet(
        counts=integer_counts, metadata=metadata, design=design, n_cpus=n_cpus, quiet=True
    )
    coeff = (
        coefficient_name(list(dds.obsm["design_matrix"].columns), variable, tested)
        if shrink
        else None
    )
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


def vst(dds: DeseqDataSet, use_design: bool = False) -> pd.DataFrame:
    """Variance-stabilized expression (samples x genes) from a pydeseq2 object."""
    dds.vst(use_design=use_design)
    return pd.DataFrame(dds.layers["vst_counts"], index=dds.obs_names, columns=dds.var_names)
