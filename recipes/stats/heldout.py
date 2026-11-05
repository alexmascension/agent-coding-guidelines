"""Recipe: held-out reference splits for empirical power and replication estimates.

Guideline: ``STATISTICS_GUIDELINES.md`` §10 (resampling, ``STAT-RES-008``): an estimate of power or
replication obtained by subsampling a large dataset defines its reference result on units disjoint
from those subsampled, with the split stratified by the variable of interest and repeated over
several random splits.

``heldout_splits`` returns, for each split, the reference units and the validation units. Within
each stratum, ``round(reference_fraction * n)`` units go to the reference (at least one, and at
least one left for validation). Seeds of the splits derive from ``seed`` through
``numpy.random.SeedSequence``, so each split is reproducible on its own.
"""

import numpy as np

REFERENCE_FRACTION = 0.5  # share of units in the reference; the split ratio is reported


def heldout_splits(
    units: np.ndarray,
    strata: np.ndarray | None = None,
    n_splits: int = 10,
    reference_fraction: float = REFERENCE_FRACTION,
    seed: int = 0,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Split experimental units into disjoint reference and validation sets, repeatedly.

    Parameters
    ----------
    units : np.ndarray
        Unit identifiers (donors), one per unit, unique.
    strata : np.ndarray | None
        Stratum of each unit (the compared groups, e.g. sex or condition); None for one stratum.
    n_splits : int
        Number of random splits.
    reference_fraction : float
        Share of units of each stratum that go to the reference, in (0, 1).
    seed : int
        Seed from which the seed of every split is derived.

    Returns
    -------
    list[tuple[np.ndarray, np.ndarray]]
        ``(reference_units, validation_units)`` per split; together they contain every unit once.

    Raises
    ------
    ValueError
        If units are repeated, lengths differ, strata have missing values, ``reference_fraction``
        is outside (0, 1), ``n_splits`` is below 1 or a stratum has fewer than 2 units.
    """
    units = np.asarray(units)
    if np.unique(units).size != units.size:
        raise ValueError("units must be unique (one entry per experimental unit)")
    if not 0 < reference_fraction < 1:
        raise ValueError(f"reference_fraction must be in (0, 1), got {reference_fraction}")
    if n_splits < 1:
        raise ValueError(f"n_splits must be at least 1, got {n_splits}")
    strata = (
        np.zeros(units.size, dtype=int) if strata is None else np.asarray(strata, dtype=object)
    )
    if strata.shape != units.shape:
        raise ValueError(f"units and strata differ in length: {units.shape} vs {strata.shape}")
    if any(s is None or (isinstance(s, float) and np.isnan(s)) for s in strata):
        raise ValueError("strata contain missing values; assign them explicitly")
    levels = sorted(set(strata.tolist()), key=str)
    members = {level: np.flatnonzero(strata == level) for level in levels}
    small = {str(level): idx.size for level, idx in members.items() if idx.size < 2}
    if small:
        raise ValueError(f"strata with fewer than 2 units: {small}")
    splits = []
    for child in np.random.SeedSequence(seed).spawn(n_splits):
        rng = np.random.default_rng(child)
        in_reference = np.zeros(units.size, dtype=bool)
        for idx in members.values():
            k = int(np.floor(reference_fraction * idx.size + 0.5))
            k = min(max(k, 1), idx.size - 1)
            in_reference[rng.choice(idx, size=k, replace=False)] = True
        splits.append((units[in_reference], units[~in_reference]))
    return splits
