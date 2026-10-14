"""Recipe: detect gene-length bias, per sample and in DE results.

Guideline: ``RNASEQ_BULK_GUIDELINES.md`` §5 (filtering and normalization) [Mandelboum 2019].

Procedure
---------
* ``sample_length_slopes``: for each sample, regress the log ratio of its normalized expression to the gene-wise
  geometric mean across samples on log gene length. A slope that differs between samples is a sample-specific length
  bias; median-of-ratios, TMM and RPKM do not remove it.
* ``slope_association``: whether slopes differ by condition (Mann-Whitney / Kruskal-Wallis). If they do, length bias
  is confounded with the condition and will produce false enrichment of short-gene (ribosomal) and long-gene
  (extracellular matrix) sets.
* ``lfc_length_trend``: Spearman correlation between log2 fold change and log gene length, with the median log2 fold
  change per length decile (the table behind the plot of the guide).

Defaults: ``MIN_MEAN_COUNT = 10`` (genes used for slopes), ``N_LENGTH_BINS = 10``. Gene length is the effective length
used for quantification (tximport ``length`` matrix averaged across samples, or the nf-core length table).
"""

import numpy as np
import pandas as pd
import polars as pl
from scipy import stats

MIN_MEAN_COUNT = 10
N_LENGTH_BINS = 10


def sample_length_slopes(
    normalized_counts: pd.DataFrame, gene_length: pd.Series, min_mean: float = MIN_MEAN_COUNT
) -> pl.DataFrame:
    """Slope of log expression ratio on log length, per sample.

    Parameters
    ----------
    normalized_counts : pd.DataFrame
        Samples x genes, counts divided by size factors.
    gene_length : pd.Series
        Gene length (bp) indexed by gene ID.
    min_mean : float
        Minimum mean normalized count of a gene to be used.

    Returns
    -------
    pl.DataFrame
        ``sample``, ``slope``, ``slope_se``, ``n_genes``.
    """
    genes = normalized_counts.columns.intersection(gene_length.index)
    x = normalized_counts[genes]
    keep = (x.mean(axis=0) >= min_mean) & (x > 0).all(axis=0)
    x = x.loc[:, keep]
    log_x = np.log2(x.to_numpy(dtype=float))
    ratio = log_x - log_x.mean(axis=0, keepdims=True)
    log_length = np.log2(gene_length[x.columns].to_numpy(dtype=float))
    rows = []
    for i, sample in enumerate(x.index):
        fit = stats.linregress(log_length, ratio[i])
        rows.append(
            {
                "sample": str(sample),
                "slope": float(fit.slope),
                "slope_se": float(fit.stderr),
                "n_genes": int(x.shape[1]),
            }
        )
    return pl.DataFrame(rows)


def slope_association(slopes: pl.DataFrame, metadata: pd.DataFrame, condition_col: str) -> dict[str, float]:
    """Test whether per-sample length slopes differ between conditions (rank-based)."""
    frame = slopes.to_pandas().set_index("sample").join(metadata[[condition_col]].astype(str), how="inner")
    groups = [g["slope"].to_numpy() for _, g in frame.groupby(condition_col)]
    if len(groups) == 2:
        result = stats.mannwhitneyu(*groups, alternative="two-sided")
        test = "mann-whitney"
    else:
        result = stats.kruskal(*groups)
        test = "kruskal-wallis"
    return {"test": test, "statistic": float(result.statistic), "p_value": float(result.pvalue)}


def lfc_length_trend(
    results: pl.DataFrame, gene_length: pd.Series, n_bins: int = N_LENGTH_BINS, lfc_col: str = "log2_fold_change"
) -> tuple[dict[str, float], pl.DataFrame]:
    """Spearman correlation of log2 fold change with log length, and median log2 fold change per length bin."""
    frame = results.select("gene_id", lfc_col).drop_nulls().to_pandas().set_index("gene_id")
    frame = frame.join(gene_length.rename("length"), how="inner").dropna()
    rho, p = stats.spearmanr(np.log2(frame["length"]), frame[lfc_col])
    frame["length_bin"] = pd.qcut(np.log2(frame["length"]), q=n_bins, labels=False, duplicates="drop")
    binned = (
        frame.groupby("length_bin")
        .agg(median_length=("length", "median"), median_lfc=(lfc_col, "median"), n_genes=(lfc_col, "size"))
        .reset_index()
    )
    return {"spearman_rho": float(rho), "p_value": float(p)}, pl.from_pandas(binned)
