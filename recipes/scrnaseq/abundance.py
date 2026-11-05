"""Recipe: cell counts per sample and cell type with explicit zeros (input of abundance methods).

Guideline: ``SCRNASEQ_GUIDELINES.md`` §11 (differential abundance). This prepares inputs; the tests
themselves (propeller, scCODA, sccomp, Milo) are run with their packages.

``cell_count_table`` returns a long table with every (sample, cell type) combination, including
zeros, the sample total and the proportion. A cell type absent from a sample is information (count
0), not a missing value.

Samples or cell types that lost all their cells (e.g. after QC) cannot be reconstructed from
``obs``: pass the sample manifest and the list of cell types of the study. Samples with no cells
get ``sample_status = "no_cells"`` and a null proportion; they MUST NOT enter a compositional
analysis, and their exclusion is recorded.

Contracts checked before counting:

* the denominator is every cell of the sample: a ``cell_types`` list that omits a type observed in
  ``obs`` raises, unless ``subcomposition=True`` declares that the proportions are within the
  listed types (the ``denominator`` column says which);
* sample identifiers of ``obs`` and of the manifest index are validated before they are compared
  as strings (``recipes/scrnaseq/units.py``);
* the manifest does not overwrite ``obs``: a covariate whose value differs between the two for a
  sample raises (correct the source that is wrong, then call again).

This module does not import a DE package (``pydeseq2`` is not needed).
"""

import pandas as pd
import polars as pl

from recipes.scrnaseq.units import conflicting_values, unit_ids, unit_metadata

DENOMINATORS = {False: "all_cell_types", True: "listed_cell_types"}


def cell_count_table(
    obs: pd.DataFrame,
    sample_key: str,
    celltype_key: str,
    covariates: list[str],
    manifest: pd.DataFrame | None = None,
    cell_types: list[str] | None = None,
    subcomposition: bool = False,
) -> pl.DataFrame:
    """Long table of counts per sample and cell type, with explicit zeros.

    Parameters
    ----------
    obs : pd.DataFrame
        ``adata.obs`` after QC and doublet removal.
    sample_key, celltype_key : str
        Columns with the sample (replicate unit) and the annotation.
    covariates : list[str]
        Sample-level columns to carry (condition, batch, age, ...); dtypes are preserved.
    manifest : pd.DataFrame | None
        All samples of the study (index = sample ID, columns include ``covariates``); None uses
        samples in ``obs``. Its values MUST agree with ``obs`` for the samples in both.
    cell_types : list[str] | None
        All cell types of the study (types that lost all their cells included); None uses those in
        ``obs``. It MUST contain every type observed in ``obs`` unless ``subcomposition`` is True.
    subcomposition : bool
        Declare that ``cell_types`` is a subset and that proportions are within it: cells of other
        types are left out of ``n_cells_sample``, and ``denominator`` is ``listed_cell_types``.

    Returns
    -------
    pl.DataFrame
        ``sample``, ``cell_type``, ``n_cells``, ``n_cells_sample`` (the denominator),
        ``proportion`` (null for samples without cells of the counted types), ``sample_status``
        (``ok`` or ``no_cells``), ``denominator`` (``all_cell_types`` or ``listed_cell_types``),
        covariates.

    Raises
    ------
    ValueError
        If sample identifiers are missing or collide as strings, covariates have missing values or
        vary within a sample, ``obs`` has samples absent from the manifest, the manifest index is
        invalid, a covariate differs between ``obs`` and the manifest, ``cell_types`` omits
        observed types without ``subcomposition``, or ``subcomposition`` is set without
        ``cell_types``.
    """
    meta = unit_metadata(obs, sample_key, covariates)
    if manifest is not None:
        index = manifest.index.to_series(index=range(len(manifest)), name="sample")
        manifest_ids = unit_ids(index.to_frame(), "sample")
        if manifest_ids.duplicated().any():
            repeated = sorted(set(manifest_ids[manifest_ids.duplicated()]))
            raise ValueError(f"samples repeated in the manifest index: {repeated[:5]}")
        manifest = manifest.copy()
        manifest.index = pd.Index(manifest_ids.to_numpy())
        unknown = sorted(set(meta.index) - set(manifest.index))
        if unknown:
            raise ValueError(f"samples in obs but not in the manifest: {unknown[:5]}")
        absent = [c for c in covariates if c not in manifest.columns]
        if absent:
            raise ValueError(f"covariates not in the manifest: {absent}")
        conflicts = conflicting_values(meta, manifest, covariates)
        if conflicts:
            shown = {column: units[:5] for column, units in conflicts.items()}
            raise ValueError(
                f"manifest and obs disagree on {shown} (column: samples); correct the wrong "
                "source before counting: the manifest does not overwrite obs"
            )
        meta = manifest[covariates]
    if subcomposition and cell_types is None:
        raise ValueError("subcomposition=True needs the list of cell_types it is defined on")
    observed_types = set(obs[celltype_key].astype(str))
    if cell_types is not None and not subcomposition:
        omitted = sorted(observed_types - {str(t) for t in cell_types})
        if omitted:
            raise ValueError(
                f"cell_types omits observed types {omitted[:5]}: proportions would not be of the "
                "sample total; add them, or pass subcomposition=True to count within the listed "
                "types"
            )
    types = sorted({str(t) for t in cell_types} if cell_types is not None else observed_types)
    observed = (
        obs.assign(sample=obs[sample_key].astype(str), cell_type=obs[celltype_key].astype(str))
        .groupby(["sample", "cell_type"], observed=True)
        .size()
    )
    grid = pd.MultiIndex.from_product([meta.index, types], names=["sample", "cell_type"])
    counts = observed.reindex(grid, fill_value=0).rename("n_cells").reset_index()
    counts["n_cells_sample"] = counts.groupby("sample")["n_cells"].transform("sum")
    counts["proportion"] = counts["n_cells"] / counts["n_cells_sample"].where(
        counts["n_cells_sample"] > 0
    )
    counts["sample_status"] = (counts["n_cells_sample"] > 0).map({True: "ok", False: "no_cells"})
    counts["denominator"] = DENOMINATORS[bool(subcomposition)]
    table = counts.merge(meta, left_on="sample", right_index=True, how="left")
    return pl.from_pandas(table).sort("sample", "cell_type")


def wide_counts(table: pl.DataFrame) -> pd.DataFrame:
    """Return the samples x cell types count matrix (samples with ``sample_status == "ok"``)."""
    usable = table.filter(pl.col("sample_status") == "ok")
    wide = usable.pivot(on="cell_type", index="sample", values="n_cells").sort("sample")
    return wide.to_pandas().set_index("sample")
