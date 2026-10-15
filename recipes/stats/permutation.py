"""Recipe: permutation of exchangeable units within strata, and permutation p-values.

Guideline: ``STATISTICS_GUIDELINES.md`` §10 (resampling): permute the labels of experimental units
(donors, samples),
within strata (batch, pair) when the design has them; $p = (b + 1)/(m + 1)$.

The caller is responsible for exchangeability: permuting labels within strata preserves the stratum
structure only.
It does not preserve the association of the labels with other covariates (age, sex, RNA quality);
when those are
associated with the labels, the permutation null is not the null of the adjusted model.

Two-sided p-values follow the SciPy convention (``scipy.stats.permutation_test``): twice the
smaller one-sided
p-value, capped at 1. Comparing absolute values is valid only for statistics centred at 0 under the
null, so it is
not used.

``scipy.stats.permutation_test`` covers unstratified tests; use this module when labels must be
permuted within strata
or when the statistic is an expensive pipeline (e.g. a full DE run,
``recipes/bulk/fdr_calibration.py``).
"""

from collections.abc import Callable

import numpy as np

ALTERNATIVES = ("two-sided", "greater", "less")


def permute_within_strata(
    labels: np.ndarray, strata: np.ndarray | None, rng: np.random.Generator
) -> np.ndarray:
    """Return a permutation of ``labels`` that shuffles only within each stratum.

    Parameters
    ----------
    labels : np.ndarray
        One label per experimental unit (not per cell).
    strata : np.ndarray | None
        Stratum per unit (batch, pair); None permutes freely. Missing strata are not allowed.
    rng : np.random.Generator
        Random generator (seeded by the caller).

    Returns
    -------
    np.ndarray
        Permuted labels; the label counts within every stratum are unchanged.

    Raises
    ------
    ValueError
        If lengths differ or strata contain missing values.
    """
    labels = np.asarray(labels)
    if strata is None:
        return rng.permutation(labels)
    strata = np.asarray(strata, dtype=object)
    if strata.shape != labels.shape:
        raise ValueError(f"labels and strata differ in length: {labels.shape} vs {strata.shape}")
    if any(s is None or (isinstance(s, float) and np.isnan(s)) for s in strata):
        raise ValueError("strata contain missing values; assign them explicitly")
    permuted = labels.copy()
    for stratum in set(strata.tolist()):
        idx = np.flatnonzero(strata == stratum)
        permuted[idx] = labels[rng.permutation(idx)]
    return permuted


def _finite_scalar(value: object, what: str) -> float:
    result = float(np.asarray(value, dtype=float))
    if not np.isfinite(result):
        raise ValueError(
            f"{what} is not finite ({result}); the statistic must return a finite scalar"
        )
    return result


def permutation_pvalue(
    statistic: Callable[[np.ndarray], float],
    labels: np.ndarray,
    strata: np.ndarray | None = None,
    n_permutations: int = 999,
    seed: int = 0,
    alternative: str = "two-sided",
) -> tuple[float, float, np.ndarray]:
    r"""Permutation p-value $p = (b + 1)/(m + 1)$ for a statistic of the unit labels.

    Parameters
    ----------
    statistic : Callable[[np.ndarray], float]
        Function of the label vector returning a finite scalar (data are closed over by the
        caller).
    labels : np.ndarray
        Observed labels of the experimental units.
    strata : np.ndarray | None
        Strata within which labels are exchangeable.
    n_permutations : int
        Number of random permutations $m$ (at least 1); the smallest attainable p-value is
        $1/(m+1)$
        (two-sided: $2/(m+1)$).
    seed : int
        Seed of the random generator.
    alternative : {"two-sided", "greater", "less"}
        ``greater``: large statistics are extreme; ``two-sided``: $\min(1, 2\min(p_{greater},
        p_{less}))$.

    Returns
    -------
    p_value : float
    observed : float
    null : np.ndarray
        Statistic under each permutation.

    Raises
    ------
    ValueError
        If the statistic is not finite (observed or any permutation), ``alternative`` is unknown,
        or
        ``n_permutations < 1``.
    """
    if alternative not in ALTERNATIVES:
        raise ValueError(f"unknown alternative {alternative!r}; use one of {ALTERNATIVES}")
    if int(n_permutations) < 1:
        raise ValueError(f"n_permutations must be >= 1, got {n_permutations}")
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)
    observed = _finite_scalar(statistic(labels), "observed statistic")
    null = np.array(
        [
            _finite_scalar(
                statistic(permute_within_strata(labels, strata, rng)),
                f"statistic of permutation {i}",
            )
            for i in range(n_permutations)
        ]
    )
    p_greater = (np.sum(null >= observed) + 1) / (n_permutations + 1)
    p_less = (np.sum(null <= observed) + 1) / (n_permutations + 1)
    if alternative == "greater":
        p = p_greater
    elif alternative == "less":
        p = p_less
    else:
        p = min(1.0, 2 * min(p_greater, p_less))
    return float(p), observed, null


def binomial_upper_bound(n_events: int, n_trials: int, level: float = 0.95) -> float:
    """One-sided Clopper-Pearson upper bound for a probability from ``n_events``/``n_trials``.

    With 0 events in 20 trials the 95% bound is about 0.139: a small number of permutations cannot
    show that a
    probability is below 0.05.
    """
    from scipy.stats import beta

    if not 0 <= n_events <= n_trials or n_trials < 1:
        raise ValueError(f"invalid counts: {n_events} events in {n_trials} trials")
    if n_events == n_trials:
        return 1.0
    return float(beta.ppf(level, n_events + 1, n_trials - n_events))
