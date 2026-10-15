"""Recipe: sensitivity of DE to the number k of hidden factors; residual-PCA estimator.

Guideline: ``RNASEQ_BULK_GUIDELINES.md`` §7 (known and hidden factors).

Procedure
---------
* ``k_sensitivity`` (the main tool): run the DE pipeline for each k in a grid, adding the first k
  factors as
  covariates, and report the number of discoveries, the overlap (Jaccard) of the discovery set with
  the one at the
  reference k, and the correlation of log2 fold changes. Factors may come from any method (qSVA,
  RUV, SVA computed in
  R and read back, or ``residual_pcs``).
* ``residual_pcs``: principal components of the expression matrix after regressing out the full
  known design.

What ``residual_pcs`` is and is not
-----------------------------------
It is an exploratory estimator for sensitivity analyses, not an implementation of qSVA, RUV or SVA.
By construction
its factors are orthogonal to the known design, so they cannot absorb the condition. The same
construction means
that the part of a technical factor aligned with the condition stays in the condition effect:
residual PCs do not
correct confounding between a hidden technical factor and the condition. The factors are also
estimated from the
genes that are then tested, and that uncertainty is not propagated. When the guide recommends qSVA
(postmortem
brain), RUV (controls or replicates) or SVA, use those methods and pass their factors to
``k_sensitivity``.
"""

import logging
from collections.abc import Callable

import numpy as np
import pandas as pd
import polars as pl

logger = logging.getLogger(__name__)

ALPHA = 0.05
SINGULAR_TOLERANCE = (
    1e-8  # relative to the largest singular value; smaller components are numerical noise
)


def residual_pcs(
    vst: pd.DataFrame, design: pd.DataFrame, k_max: int, n_top_genes: int | None = 2000
) -> pd.DataFrame:
    """Principal components of expression residuals after regressing out the known design.

    Parameters
    ----------
    vst : pd.DataFrame
        Samples x genes, variance-stabilized expression.
    design : pd.DataFrame
        Samples x covariates of the known model (condition included); categorical columns are
        dummy-coded.
    k_max : int
        Number of factors returned (>= 1).
    n_top_genes : int | None
        Use the most variable residual genes (None uses all).

    Returns
    -------
    pd.DataFrame
        Samples x ``hf1..hf{k_max}``, standardized.

    Raises
    ------
    ValueError
        If indices differ, ``k_max`` is not positive, or ``k_max`` exceeds the number of
        non-degenerate residual
        components (limited by samples minus the rank of the design, and by the number of genes).
    """
    if not vst.index.equals(design.index):
        raise ValueError("vst and design must have the same sample index, in the same order")
    if k_max < 1:
        raise ValueError(f"k_max must be >= 1, got {k_max}")
    x_design = pd.get_dummies(design, drop_first=True, dtype=float)
    x_design.insert(0, "intercept", 1.0)
    x = x_design.to_numpy()
    design_rank = int(np.linalg.matrix_rank(x))
    y = vst.to_numpy(dtype=float)
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    residuals = y - x @ beta
    if n_top_genes is not None and n_top_genes < residuals.shape[1]:
        top = np.argsort(residuals.var(axis=0))[::-1][:n_top_genes]
        residuals = residuals[:, top]
    u, s, _ = np.linalg.svd(residuals - residuals.mean(axis=0), full_matrices=False)
    usable = int(np.sum(s > SINGULAR_TOLERANCE * s[0])) if s.size and s[0] > 0 else 0
    limit = min(usable, vst.shape[0] - design_rank, residuals.shape[1])
    if k_max > limit:
        raise ValueError(
            f"k_max={k_max} exceeds the available residual components ({limit}: "
            f"{vst.shape[0]} samples, design rank {design_rank}, {residuals.shape[1]} genes, "
            f"{usable} non-degenerate components)"
        )
    factors = u[:, :k_max] * s[:k_max]
    factors = (factors - factors.mean(axis=0)) / factors.std(axis=0)
    return pd.DataFrame(factors, index=vst.index, columns=[f"hf{i + 1}" for i in range(k_max)])


def k_sensitivity(
    de_pipeline: Callable[[pd.DataFrame, list[str]], pl.DataFrame],
    metadata: pd.DataFrame,
    factors: pd.DataFrame,
    k_grid: list[int],
    k_reference: int,
    alpha: float = ALPHA,
) -> pl.DataFrame:
    """DE result as a function of the number of hidden factors.

    Parameters
    ----------
    de_pipeline : Callable[[pd.DataFrame, list[str]], pl.DataFrame]
        Takes metadata (with factor columns added) and the list of factor column names to add to
        the design; returns
        a table with ``gene_id``, ``log2_fold_change`` and ``p_adjusted``.
    metadata : pd.DataFrame
        Sample metadata.
    factors : pd.DataFrame
        Samples x factors (``residual_pcs`` or an external method); every metadata sample must be
        present.
    k_grid : list[int]
        Numbers of factors to test (non-negative; include 0 and the method's estimate).
    k_reference : int
        k of the main analysis; overlaps are computed against it.
    alpha : float
        FDR threshold for discoveries.

    Returns
    -------
    pl.DataFrame
        ``k``, ``n_discoveries``, ``jaccard_vs_reference``, ``lfc_correlation_vs_reference``.

    Raises
    ------
    ValueError
        If ``k_reference`` is not in the grid, a k is negative or exceeds the available factors, or
        samples are
        missing from ``factors``.
    """
    if k_reference not in k_grid:
        raise ValueError("k_reference must be in k_grid")
    if min(k_grid) < 0 or max(k_grid) > factors.shape[1]:
        raise ValueError(f"k must be in [0, {factors.shape[1]}], got {sorted(k_grid)}")
    missing = metadata.index.difference(factors.index)
    if len(missing):
        raise ValueError(f"samples without factor values: {list(missing[:5])}")
    factors = factors.loc[metadata.index]
    results = {}
    for k in sorted(set(k_grid)):
        columns = list(factors.columns[:k])
        meta = metadata.join(factors[columns]) if k else metadata.copy()
        results[k] = de_pipeline(meta, columns)
        logger.info(f"k={k}: {int((results[k]['p_adjusted'] < alpha).sum())} discoveries")
    reference = results[k_reference]
    ref_set = set(reference.filter(pl.col("p_adjusted") < alpha)["gene_id"].to_list())
    rows = []
    for k, table in results.items():
        found = set(table.filter(pl.col("p_adjusted") < alpha)["gene_id"].to_list())
        union = found | ref_set
        joined = (
            table.select("gene_id", "log2_fold_change")
            .join(
                reference.select("gene_id", pl.col("log2_fold_change").alias("ref_lfc")),
                on="gene_id",
            )
            .drop_nulls()
        )
        lfc, ref = joined["log2_fold_change"].to_numpy(), joined["ref_lfc"].to_numpy()
        correlation = (
            float(np.corrcoef(lfc, ref)[0, 1]) if lfc.std() > 0 and ref.std() > 0 else float("nan")
        )
        rows.append(
            {
                "k": k,
                "n_discoveries": len(found),
                "jaccard_vs_reference": len(found & ref_set) / len(union) if union else 1.0,
                "lfc_correlation_vs_reference": correlation,
            }
        )
    return pl.DataFrame(rows)
