"""Recipe: hidden factors with the condition protected, and sensitivity of the DE result to their number k.

Guideline: ``RNASEQ_BULK_GUIDELINES.md`` §7 (known and hidden factors).

Procedure
---------
* ``residual_pcs``: principal components of the VST matrix after regressing out the full known design (condition
  and known covariates). Because the condition is in the regression, the factors cannot absorb it; this is the
  protection that plain PCA of expression lacks. It is a Python estimator in the spirit of RUVr/svaseq; when the guide
  prefers qSVA, RUV or SVA, compute the factors with those methods (R script) and pass them to ``k_sensitivity``.
* ``k_sensitivity``: run the DE pipeline for each k in a grid, adding the first k factors as covariates, and report
  the number of discoveries, the overlap (Jaccard) of the discovery set with the one at the reference k, and the
  correlation of log2 fold changes. Conclusions that change across neighbouring k depend on the choice.

Factors are continuous covariates; standardize them (done here) so that pydeseq2 fits are well conditioned.
"""

import logging
from collections.abc import Callable

import numpy as np
import pandas as pd
import polars as pl

logger = logging.getLogger(__name__)

ALPHA = 0.05


def residual_pcs(vst: pd.DataFrame, design: pd.DataFrame, k_max: int, n_top_genes: int | None = 2000) -> pd.DataFrame:
    """Principal components of VST residuals after regressing out the known design.

    Parameters
    ----------
    vst : pd.DataFrame
        Samples x genes, variance-stabilized expression.
    design : pd.DataFrame
        Samples x covariates of the known model (condition included); categorical columns are dummy-coded.
    k_max : int
        Number of factors returned.
    n_top_genes : int | None
        Use the most variable residual genes (None uses all).

    Returns
    -------
    pd.DataFrame
        Samples x ``hf1..hf{k_max}``, standardized.
    """
    if not vst.index.equals(design.index):
        raise ValueError("vst and design must have the same sample index, in the same order")
    x_design = pd.get_dummies(design, drop_first=True, dtype=float)
    x_design.insert(0, "intercept", 1.0)
    y = vst.to_numpy(dtype=float)
    beta, *_ = np.linalg.lstsq(x_design.to_numpy(), y, rcond=None)
    residuals = y - x_design.to_numpy() @ beta
    if n_top_genes is not None and n_top_genes < residuals.shape[1]:
        top = np.argsort(residuals.var(axis=0))[::-1][:n_top_genes]
        residuals = residuals[:, top]
    residuals -= residuals.mean(axis=0)
    u, s, _ = np.linalg.svd(residuals, full_matrices=False)
    max_rank = residuals.shape[0] - x_design.shape[1]
    if k_max > max_rank:
        raise ValueError(f"k_max={k_max} exceeds the residual degrees of freedom ({max_rank})")
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
        Takes metadata (with factor columns added) and the list of factor column names to add to the design; returns
        a table with ``gene_id``, ``log2_fold_change`` and ``p_adjusted``.
    metadata : pd.DataFrame
        Sample metadata.
    factors : pd.DataFrame
        Samples x factors (``residual_pcs`` or an external method), same index as ``metadata``.
    k_grid : list[int]
        Numbers of factors to test (include 0 and the method's estimate).
    k_reference : int
        k of the main analysis; overlaps are computed against it.
    alpha : float
        FDR threshold for discoveries.

    Returns
    -------
    pl.DataFrame
        ``k``, ``n_discoveries``, ``jaccard_vs_reference``, ``lfc_correlation_vs_reference``.
    """
    if k_reference not in k_grid:
        raise ValueError("k_reference must be in k_grid")
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
            .join(reference.select("gene_id", pl.col("log2_fold_change").alias("ref_lfc")), on="gene_id")
            .drop_nulls()
        )
        rows.append(
            {
                "k": k,
                "n_discoveries": len(found),
                "jaccard_vs_reference": len(found & ref_set) / len(union) if union else 1.0,
                "lfc_correlation_vs_reference": float(np.corrcoef(joined["log2_fold_change"], joined["ref_lfc"])[0, 1]),
            }
        )
    return pl.DataFrame(rows)
