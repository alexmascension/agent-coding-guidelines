r"""Recipe: check FDR calibration of a DE pipeline by re-running it with permuted condition labels.

Guideline: ``RNASEQ_BULK_GUIDELINES.md`` §8 (differential expression; large human cohorts) [Li 2022].

Procedure
---------
1. Permute the condition labels of samples (``recipes.stats.permutation.permute_within_strata``), within strata such
   as batch or sequencing platform so that the permuted design keeps the same confounding structure.
2. Re-run the complete DE pipeline (same filter, design and covariates) on each permuted design.
3. Under permuted labels no gene is truly DE, so every discovery is false. With FDR control at level $\alpha$, the
   fraction of permutations with at least one discovery should be at most about $\alpha$, and most permutations should
   report zero genes. Many discoveries in many permutations mean the pipeline is anti-conservative for this cohort
   (e.g. outliers, unmodelled heterogeneity), and its real result is not trustworthy as is.

The DE pipeline is passed as a function so the same check applies to pydeseq2, edgeR/limma scripts or pseudobulk.
Each permutation costs one full DE run; ``N_PERMUTATIONS = 20`` is a minimum for a qualitative check.
"""

import logging
from collections.abc import Callable

import numpy as np
import pandas as pd
import polars as pl

from recipes.stats.permutation import permute_within_strata

logger = logging.getLogger(__name__)

N_PERMUTATIONS = 20
ALPHA = 0.05


def fdr_calibration(
    de_pipeline: Callable[[pd.DataFrame], np.ndarray],
    metadata: pd.DataFrame,
    condition_col: str,
    strata_col: str | None = None,
    n_permutations: int = N_PERMUTATIONS,
    alpha: float = ALPHA,
    seed: int = 0,
) -> tuple[pl.DataFrame, dict[str, float]]:
    """Count discoveries of a DE pipeline under permuted condition labels.

    Parameters
    ----------
    de_pipeline : Callable[[pd.DataFrame], np.ndarray]
        Takes the sample metadata (with possibly permuted ``condition_col``) and returns adjusted p-values for all
        tested genes. It must apply the same filtering and design as the real analysis.
    metadata : pd.DataFrame
        Sample metadata of the real analysis.
    condition_col : str
        Column with the condition labels to permute.
    strata_col : str | None
        Column within which labels are permuted (batch, platform); None permutes freely.
    n_permutations : int
        Number of permuted runs.
    alpha : float
        FDR level of the real analysis.
    seed : int
        Seed of the permutations.

    Returns
    -------
    per_permutation : pl.DataFrame
        ``permutation``, ``n_discoveries``, ``n_tested``, ``fraction_label_changed`` (share of samples whose label
        changed; permutations that barely change labels are not informative).
    summary : dict[str, float]
        ``fraction_with_discoveries``, ``median_discoveries``, ``max_discoveries``, ``alpha``.
    """
    rng = np.random.default_rng(seed)
    labels = metadata[condition_col].to_numpy()
    strata = metadata[strata_col].to_numpy() if strata_col else None
    rows = []
    for i in range(n_permutations):
        permuted = metadata.copy()
        new_labels = permute_within_strata(labels, strata, rng)
        permuted[condition_col] = new_labels
        p_adjusted = np.asarray(de_pipeline(permuted), dtype=float)
        rows.append(
            {
                "permutation": i,
                "n_discoveries": int(np.sum(np.nan_to_num(p_adjusted, nan=1.0) < alpha)),
                "n_tested": int(np.sum(np.isfinite(p_adjusted))),
                "fraction_label_changed": float(np.mean(new_labels != labels)),
            }
        )
        logger.info(f"permutation {i + 1}/{n_permutations}: {rows[-1]['n_discoveries']} discoveries")
    table = pl.DataFrame(rows)
    summary = {
        "fraction_with_discoveries": float((table["n_discoveries"] > 0).mean()),
        "median_discoveries": float(table["n_discoveries"].median()),
        "max_discoveries": float(table["n_discoveries"].max()),
        "alpha": alpha,
    }
    return table, summary
