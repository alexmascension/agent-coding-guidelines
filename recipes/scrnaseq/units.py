"""Recipe: identifiers and metadata of experimental units (samples, donors), validated.

Guidelines: ``SCRNASEQ_GUIDELINES.md`` §2 (design) and §10-11; ``CODING_GUIDELINES.md``
``CODE-DATA-008`` (keys are validated before they are used).

Unit identifiers are compared as strings downstream. Before converting them, ``unit_ids`` rejects
missing identifiers (they would become one unit called ``"nan"``) and distinct values that become
the same string (``1`` and ``"1"`` would be merged into one unit). When a project decides that such
values are the same unit, it canonicalizes them before calling and records that decision.

This module needs no DE package: ``recipes/scrnaseq/abundance.py`` imports it without
``pydeseq2``.
"""

import numpy as np
import pandas as pd


def unit_ids(obs: pd.DataFrame, key: str) -> pd.Series:
    """Return ``obs[key]`` as strings after checking that the identifiers are well defined.

    Parameters
    ----------
    obs : pd.DataFrame
        Cell (or sample) metadata.
    key : str
        Column with the identifiers (unit, sample, donor).

    Returns
    -------
    pd.Series
        The identifiers as strings, with the index of ``obs``.

    Raises
    ------
    ValueError
        If the column is absent, has missing values, or distinct values become the same string.
    """
    if key not in obs.columns:
        raise ValueError(f"column {key!r} not in obs")
    values = obs[key]
    missing = values.isna().to_numpy()
    if missing.any():
        raise ValueError(
            f"{int(missing.sum())} rows with a missing {key!r}; assign the identifier or remove "
            "the rows before aggregation"
        )
    as_str = values.astype(str)
    pairs = pd.DataFrame(
        {"as_str": as_str.to_numpy(), "value": values.astype(object).to_numpy()}
    ).drop_duplicates()
    collide = pairs[pairs["as_str"].duplicated(keep=False)]
    if not collide.empty:
        examples = collide.groupby("as_str")["value"].agg(lambda v: [repr(x) for x in v])
        raise ValueError(
            f"distinct values of {key!r} become the same identifier as strings: "
            f"{examples.head(5).to_dict()}; canonicalize them before calling and record it"
        )
    return as_str


def unit_metadata(obs: pd.DataFrame, unit_key: str, covariates: list[str]) -> pd.DataFrame:
    """One row per aggregation unit with its covariates, keeping dtypes.

    Raises
    ------
    ValueError
        If the unit identifiers are missing or collide as strings (``unit_ids``), or a covariate
        has missing values or is not constant within a unit.
    """
    units = unit_ids(obs, unit_key)
    missing = [c for c in covariates if obs[c].isna().any()]
    if missing:
        raise ValueError(
            f"covariates with missing values: {missing}; resolve them before aggregation"
        )
    varying = [
        c
        for c in covariates
        if (obs.groupby(units, observed=True)[c].nunique(dropna=False) > 1).any()
    ]
    if varying:
        raise ValueError(f"covariates not constant within {unit_key}: {varying}")
    meta = obs.assign(_unit=units).drop_duplicates("_unit").set_index("_unit")[covariates]
    meta.index.name = unit_key
    for column in covariates:
        if isinstance(meta[column].dtype, pd.CategoricalDtype):
            meta[column] = meta[column].cat.remove_unused_categories()
    return meta


def _equal_values(a: object, b: object) -> bool:
    """Whether two metadata values are the same (``30`` and ``30.0``, or ``"30"`` and ``30``)."""
    if pd.isna(a) or pd.isna(b):
        return bool(pd.isna(a) and pd.isna(b))
    try:
        if a == b:
            return True
    except (TypeError, ValueError):
        pass
    return str(a) == str(b)


def conflicting_values(
    observed: pd.DataFrame, declared: pd.DataFrame, columns: list[str]
) -> dict[str, list[str]]:
    """Units whose ``columns`` differ between two tables indexed by unit (e.g. obs and manifest).

    Returns
    -------
    dict[str, list[str]]
        Column -> units (of ``observed``) with a different value in ``declared``; empty when they
        agree.
    """
    conflicts = {}
    for column in columns:
        a = np.asarray(observed[column], dtype=object)
        b = np.asarray(declared.loc[observed.index, column], dtype=object)
        bad = [
            str(unit)
            for unit, x, y in zip(observed.index, a, b, strict=True)
            if not _equal_values(x, y)
        ]
        if bad:
            conflicts[column] = bad
    return conflicts
