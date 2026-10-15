r"""Recipe: global-null diagnostic of a DE pipeline, re-run with permuted condition labels.

Guideline: ``RNASEQ_BULK_GUIDELINES.md`` §8 (differential expression; large human cohorts) [Li
2022].

What it checks, and what it does not
------------------------------------
Re-running the full pipeline with permuted condition labels creates data in which the condition has
no effect,
**provided the permuted samples are exchangeable under the null**. Every discovery under
permutation is then false,
and with FDR control at level $\alpha$ the probability of at least one discovery should be at most
about $\alpha$
(under the global null, FDR equals that probability). Many discoveries in many permutations mean
the pipeline is
anti-conservative for this cohort (outliers, unmodelled heterogeneity).

* Exchangeability MUST be justified by the caller. Permuting within strata (batch, platform) keeps
  the strata
  structure only; it does not keep the association of the condition with other covariates (age,
  sex, RNA quality,
  postmortem interval). In observational cohorts where those covariates differ by condition, the
  permuted labels
  break that association, and the diagnostic answers a different question. Alternatives: permute
  within finer strata
  (e.g. matched sets), or simulate from the fitted null model (parametric simulation).
* Passing the diagnostic does not establish FDR control when some genes are truly DE; it only shows
  that the pipeline
  does not produce discoveries from nothing.
* With few permutations the diagnostic is qualitative: the summary reports the Clopper-Pearson
  upper bound for the
  probability of any discovery (0 of 20 permutations still allows up to about 14%).

The DE pipeline is passed as a function so the same check applies to pydeseq2, edgeR/limma scripts
or pseudobulk.
Each permutation costs one full DE run; ``N_PERMUTATIONS = 20`` is a minimum for a qualitative
check.
"""

import logging
from collections.abc import Callable

import numpy as np
import pandas as pd
import polars as pl

from recipes.stats.permutation import binomial_upper_bound, permute_within_strata

logger = logging.getLogger(__name__)

N_PERMUTATIONS = 20
ALPHA = 0.05
CONFIDENCE = 0.95


def fdr_calibration(
    de_pipeline: Callable[[pd.DataFrame], np.ndarray],
    metadata: pd.DataFrame,
    condition_col: str,
    strata_col: str | None = None,
    n_permutations: int = N_PERMUTATIONS,
    alpha: float = ALPHA,
    seed: int = 0,
) -> tuple[pl.DataFrame, dict[str, float]]:
    """Count discoveries of a DE pipeline under permuted condition labels (global-null diagnostic).

    Parameters
    ----------
    de_pipeline : Callable[[pd.DataFrame], np.ndarray]
        Takes the sample metadata (with possibly permuted ``condition_col``) and returns adjusted
        p-values for all
        genes (NaN for genes not tested). It must apply the same filtering and design as the real
        analysis.
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
        ``permutation``, ``n_discoveries``, ``n_tested``, ``n_not_tested``,
        ``fraction_label_changed`` (permutations
        that barely change labels are not informative).
    summary : dict[str, float]
        ``n_permutations``, ``n_distinct_permutations``, ``fraction_with_discoveries``,
        ``upper_bound_any_discovery`` (one-sided 95% Clopper-Pearson), ``median_discoveries``,
        ``max_discoveries``,
        ``alpha``, ``strata``.

    Raises
    ------
    ValueError
        If the pipeline returns adjusted p-values outside [0, 1].
    """
    rng = np.random.default_rng(seed)
    labels = metadata[condition_col].to_numpy()
    strata = metadata[strata_col].to_numpy() if strata_col else None
    rows, seen = [], set()
    for i in range(n_permutations):
        permuted = metadata.copy()
        new_labels = permute_within_strata(labels, strata, rng)
        seen.add(tuple(new_labels.tolist()))
        permuted[condition_col] = new_labels
        p_adjusted = np.asarray(de_pipeline(permuted), dtype=float)
        finite = np.isfinite(p_adjusted)
        if np.any((p_adjusted[finite] < 0) | (p_adjusted[finite] > 1)):
            raise ValueError(f"permutation {i}: adjusted p-values outside [0, 1]")
        rows.append(
            {
                "permutation": i,
                "n_discoveries": int(np.sum(p_adjusted[finite] < alpha)),
                "n_tested": int(finite.sum()),
                "n_not_tested": int((~finite).sum()),
                "fraction_label_changed": float(np.mean(new_labels != labels)),
            }
        )
        logger.info(
            f"permutation {i + 1}/{n_permutations}: {rows[-1]['n_discoveries']} discoveries"
        )
    table = pl.DataFrame(rows)
    n_with = int((table["n_discoveries"] > 0).sum())
    summary = {
        "n_permutations": n_permutations,
        "n_distinct_permutations": len(seen),
        "fraction_with_discoveries": n_with / n_permutations,
        "upper_bound_any_discovery": binomial_upper_bound(n_with, n_permutations, CONFIDENCE),
        "median_discoveries": float(table["n_discoveries"].median()),
        "max_discoveries": float(table["n_discoveries"].max()),
        "alpha": alpha,
        "strata": strata_col,
    }
    return table, summary
