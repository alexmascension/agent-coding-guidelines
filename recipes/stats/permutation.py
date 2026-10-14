"""Recipe: permutation of exchangeable units within strata, and permutation p-values.

Guideline: ``STATISTICS_GUIDELINES.md`` §10 (resampling): permute the labels of experimental units (donors, samples),
within strata (batch, pair) when the design has them; $p = (b + 1)/(m + 1)$.

``scipy.stats.permutation_test`` covers unstratified tests; use this module when labels must be permuted within strata
or when the statistic is an expensive pipeline (e.g. a full DE run, ``recipes/bulk/fdr_calibration.py``).
"""

from collections.abc import Callable

import numpy as np


def permute_within_strata(labels: np.ndarray, strata: np.ndarray | None, rng: np.random.Generator) -> np.ndarray:
    """Return a permutation of ``labels`` that shuffles only within each stratum.

    Parameters
    ----------
    labels : np.ndarray
        One label per experimental unit (not per cell).
    strata : np.ndarray | None
        Stratum per unit (batch, pair); None permutes freely.
    rng : np.random.Generator
        Random generator (seeded by the caller).

    Returns
    -------
    np.ndarray
        Permuted labels; the label counts within every stratum are unchanged.
    """
    labels = np.asarray(labels)
    if strata is None:
        return rng.permutation(labels)
    strata = np.asarray(strata)
    if strata.shape != labels.shape:
        raise ValueError(f"labels and strata differ in length: {labels.shape} vs {strata.shape}")
    permuted = labels.copy()
    for stratum in np.unique(strata):
        idx = np.flatnonzero(strata == stratum)
        permuted[idx] = labels[rng.permutation(idx)]
    return permuted


def permutation_pvalue(
    statistic: Callable[[np.ndarray], float],
    labels: np.ndarray,
    strata: np.ndarray | None = None,
    n_permutations: int = 999,
    seed: int = 0,
    alternative: str = "two-sided",
) -> tuple[float, float, np.ndarray]:
    """Permutation p-value $p = (b + 1)/(m + 1)$ for a statistic of the unit labels.

    Parameters
    ----------
    statistic : Callable[[np.ndarray], float]
        Function of the label vector returning the test statistic (data are closed over by the caller).
    labels : np.ndarray
        Observed labels of the experimental units.
    strata : np.ndarray | None
        Strata within which labels are exchangeable.
    n_permutations : int
        Number of random permutations $m$; the smallest attainable p-value is $1/(m+1)$.
    seed : int
        Seed of the random generator.
    alternative : {"two-sided", "greater", "less"}
        Direction; two-sided compares absolute values.

    Returns
    -------
    p_value : float
    observed : float
    null : np.ndarray
        Statistic under each permutation.
    """
    rng = np.random.default_rng(seed)
    observed = float(statistic(np.asarray(labels)))
    null = np.array([statistic(permute_within_strata(labels, strata, rng)) for _ in range(n_permutations)])
    if alternative == "two-sided":
        b = int(np.sum(np.abs(null) >= abs(observed)))
    elif alternative == "greater":
        b = int(np.sum(null >= observed))
    elif alternative == "less":
        b = int(np.sum(null <= observed))
    else:
        raise ValueError(f"unknown alternative {alternative!r}")
    return (b + 1) / (n_permutations + 1), observed, null
