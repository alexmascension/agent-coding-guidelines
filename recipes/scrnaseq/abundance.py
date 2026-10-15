"""Recipe: cell counts per sample and cell type with explicit zeros (input of abundance methods).

Guideline: ``SCRNASEQ_GUIDELINES.md`` §11 (differential abundance). This prepares inputs; the tests
themselves
(propeller, scCODA, sccomp, Milo) are run with their packages.

``cell_count_table`` returns a long table with every (sample, cell type) combination, including
zeros, the sample
total and the proportion. A cell type absent from a sample is information (count 0), not a missing
value.

Samples or cell types that lost all their cells (e.g. after QC) cannot be reconstructed from
``obs``: pass the sample
manifest and the list of cell types of the study. Samples with no cells get ``sample_status =
"no_cells"`` and a null
proportion; they MUST NOT enter a compositional analysis, and their exclusion is recorded.
"""

import pandas as pd
import polars as pl

from recipes.scrnaseq.pseudobulk import unit_metadata


def cell_count_table(
    obs: pd.DataFrame,
    sample_key: str,
    celltype_key: str,
    covariates: list[str],
    manifest: pd.DataFrame | None = None,
    cell_types: list[str] | None = None,
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
        samples in ``obs``.
    cell_types : list[str] | None
        All cell types of the study; None uses those in ``obs``.

    Returns
    -------
    pl.DataFrame
        ``sample``, ``cell_type``, ``n_cells``, ``n_cells_sample``, ``proportion`` (null for
        samples without cells),
        ``sample_status`` (``ok`` or ``no_cells``), covariates.

    Raises
    ------
    ValueError
        If covariates have missing values or vary within a sample, or ``obs`` has samples absent
        from the manifest.
    """
    meta = unit_metadata(obs, sample_key, covariates)
    if manifest is not None:
        manifest = manifest.copy()
        manifest.index = manifest.index.astype(str)
        unknown = sorted(set(meta.index) - set(manifest.index))
        if unknown:
            raise ValueError(f"samples in obs but not in the manifest: {unknown[:5]}")
        meta = manifest[covariates]
    types = sorted(
        set(cell_types) if cell_types is not None else set(obs[celltype_key].astype(str))
    )
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
    table = counts.merge(meta, left_on="sample", right_index=True, how="left")
    return pl.from_pandas(table).sort("sample", "cell_type")


def wide_counts(table: pl.DataFrame) -> pd.DataFrame:
    """Return the samples x cell types count matrix (samples with ``sample_status == "ok"``)."""
    usable = table.filter(pl.col("sample_status") == "ok")
    wide = usable.pivot(on="cell_type", index="sample", values="n_cells").sort("sample")
    return wide.to_pandas().set_index("sample")
