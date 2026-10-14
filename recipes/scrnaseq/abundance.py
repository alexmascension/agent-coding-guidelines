"""Recipe: cell counts per sample and cell type with explicit zeros (input of every abundance method).

Guideline: ``SCRNASEQ_GUIDELINES.md`` §11 (differential abundance).

``cell_count_table`` returns a long table with every (sample, cell type) combination, including zeros, the sample
total and the proportion. A cell type absent from a sample is information (count 0), not a missing value; dropping
those rows biases every compositional method. ``wide_counts`` gives the samples x cell types matrix used by scCODA,
propeller or sccomp.
"""

import pandas as pd
import polars as pl


def cell_count_table(obs: pd.DataFrame, sample_key: str, celltype_key: str, covariates: list[str]) -> pl.DataFrame:
    """Long table of counts per sample and cell type, with explicit zeros.

    Parameters
    ----------
    obs : pd.DataFrame
        ``adata.obs`` after QC and doublet removal.
    sample_key, celltype_key : str
        Columns with the sample (replicate unit) and the annotation.
    covariates : list[str]
        Sample-level columns to carry (condition, batch, ...); must be constant within sample.

    Returns
    -------
    pl.DataFrame
        ``sample``, ``cell_type``, ``n_cells``, ``n_cells_sample``, ``proportion``, covariates.

    Raises
    ------
    ValueError
        If a covariate is not constant within a sample.
    """
    frame = pl.from_pandas(obs[[sample_key, celltype_key, *covariates]].astype(str).reset_index(drop=True))
    frame = frame.rename({sample_key: "sample", celltype_key: "cell_type"})
    per_sample = frame.group_by("sample").agg([pl.col(c).n_unique().alias(c) for c in covariates])
    bad = [c for c in covariates if (per_sample[c] > 1).any()]
    if bad:
        raise ValueError(f"covariates not constant within samples: {bad}")
    samples = frame.select("sample", *covariates).unique("sample")
    grid = samples.select("sample").join(frame.select("cell_type").unique(), how="cross")
    counts = frame.group_by("sample", "cell_type").agg(pl.len().alias("n_cells"))
    table = grid.join(counts, on=["sample", "cell_type"], how="left").with_columns(pl.col("n_cells").fill_null(0))
    table = table.with_columns(n_cells_sample=pl.col("n_cells").sum().over("sample"))
    table = table.with_columns(proportion=pl.col("n_cells") / pl.col("n_cells_sample"))
    return table.join(samples, on="sample").sort("sample", "cell_type")


def wide_counts(table: pl.DataFrame) -> pd.DataFrame:
    """Return the samples x cell types count matrix from ``cell_count_table``."""
    wide = table.pivot(on="cell_type", index="sample", values="n_cells").sort("sample")
    return wide.to_pandas().set_index("sample")
