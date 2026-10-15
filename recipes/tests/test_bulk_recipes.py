import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import pytest
from scipy import stats
from statsmodels.stats.multitest import multipletests

from recipes.bulk.counts_io import mean_gene_length, read_gene_table, to_integer_counts
from recipes.bulk.deseq2 import filter_genes, run_deseq2, vst
from recipes.bulk.fdr_calibration import fdr_calibration
from recipes.bulk.hidden_factors import k_sensitivity, residual_pcs
from recipes.bulk.length_bias import lfc_length_trend, sample_length_slopes, slope_association
from recipes.stats.permutation import permutation_pvalue, permute_within_strata

RSCRIPT = os.environ.get("RECIPES_RSCRIPT", "Rscript")
TXIMPORT_R = Path(__file__).resolve().parents[1] / "bulk" / "tximport_lengthscaled.R"


def simulate_bulk(
    seed=0, n_genes=500, n_per_group=6, n_de=25, fold=3.0, length_slope_case=0.0, hidden=0.0
):
    """Bulk counts: control/case x 2 batches, ``n_de`` up-regulated genes, optional length bias and hidden factor."""
    rng = np.random.default_rng(seed)
    n = 2 * n_per_group
    meta = pd.DataFrame(
        {
            "condition": ["control"] * n_per_group + ["case"] * n_per_group,
            "batch": (["b1", "b2"] * n_per_group)[:n],
        },
        index=[f"S{i}" for i in range(n)],
    )
    length = pd.Series(
        np.exp(rng.uniform(np.log(500), np.log(20000), n_genes)),
        index=[f"ENSG{i:011d}" for i in range(n_genes)],
    )
    base = rng.gamma(1.5, 100, n_genes)
    hidden_factor = rng.normal(size=n)
    loadings = rng.normal(size=n_genes) * hidden
    counts = np.empty((n, n_genes))
    for i in range(n):
        mu = base * rng.lognormal(0, 0.2)
        if meta["condition"].iloc[i] == "case":
            mu[:n_de] *= fold
            mu *= (length.to_numpy() / 2000.0) ** length_slope_case
        mu *= np.exp(loadings * hidden_factor[i])
        counts[i] = rng.negative_binomial(10, 10 / (10 + mu))
    return (
        pd.DataFrame(counts, index=meta.index, columns=length.index),
        meta,
        length,
        hidden_factor,
    )


def simple_pipeline(counts):
    """Fast stand-in DE pipeline (t-test on log CPM + BH) for permutation tests."""
    log_cpm = np.log2(counts.div(counts.sum(axis=1), axis=0) * 1e6 + 1)

    def run(meta, factor_columns=()):
        groups = meta["condition"].to_numpy()
        p = stats.ttest_ind(log_cpm[groups == "case"], log_cpm[groups == "control"], axis=0).pvalue
        return multipletests(p, method="fdr_bh")[1]

    return run


# --- permutation
# ------------------------------------------------------------------------------------------------
def test_permutation_keeps_counts_within_strata():
    rng = np.random.default_rng(0)
    labels = np.array(["a", "a", "b", "b", "a", "b"])
    strata = np.array([1, 1, 1, 2, 2, 2])
    permuted = permute_within_strata(labels, strata, rng)
    for s in (1, 2):
        assert sorted(permuted[strata == s]) == sorted(labels[strata == s])


def test_permutation_pvalue_formula_and_power():
    rng = np.random.default_rng(1)
    values = np.r_[rng.normal(0, 1, 10), rng.normal(3, 1, 10)]
    labels = np.array([0] * 10 + [1] * 10)
    p, _, null = permutation_pvalue(
        lambda lab: values[lab == 1].mean() - values[lab == 0].mean(), labels, n_permutations=199
    )
    assert 1 / 200 <= p <= 0.02
    assert null.size == 199


# --- DESeq2
# -----------------------------------------------------------------------------------------------------
def test_filter_genes_uses_smallest_group():
    counts = pd.DataFrame({"g1": [10, 10, 0, 0, 0], "g2": [10, 0, 0, 0, 0]}, index=list("abcde"))
    meta = pd.DataFrame({"condition": ["x", "x", "y", "y", "y"]}, index=list("abcde"))
    assert list(filter_genes(counts, meta, "condition")) == ["g1"]


def test_run_deseq2_recovers_planted_genes_and_threshold_test():
    counts, meta, _, _ = simulate_bulk()
    genes = filter_genes(counts, meta, "condition")
    results, dds = run_deseq2(
        counts[genes], meta, "~batch + condition", ("condition", "case", "control")
    )
    found = set(results.filter(pl.col("p_adjusted") < 0.05)["gene_id"])
    planted = set(counts.columns[:25])
    assert len(found & planted) >= 20
    assert len(found - planted) / len(found) < 0.2
    assert results.columns[:3] == ["gene_id", "base_mean", "log2_fold_change"]
    assert results["log2_fold_change_shrunk"].is_not_null().all()
    thresholded, _ = run_deseq2(
        counts[genes],
        meta,
        "~batch + condition",
        ("condition", "case", "control"),
        lfc_null=2.0,
        shrink=False,
    )
    assert (thresholded["p_adjusted"] < 0.05).sum() < len(found)
    assert vst(dds).shape == (counts.shape[0], len(genes))


# --- FDR calibration
# --------------------------------------------------------------------------------------------
def test_fdr_calibration_on_null_data_has_few_discoveries():
    counts, meta, _, _ = simulate_bulk(n_de=0, seed=4)
    table, summary = fdr_calibration(
        simple_pipeline(counts), meta, "condition", strata_col="batch", n_permutations=10
    )
    assert table.height == 10
    assert summary["median_discoveries"] == 0


# --- length bias
# ------------------------------------------------------------------------------------------------
def test_length_bias_detected_when_planted():
    counts, meta, length, _ = simulate_bulk(n_de=0, length_slope_case=0.3, seed=2)
    normalized = counts.div(counts.sum(axis=1) / counts.sum(axis=1).mean(), axis=0)
    slopes = sample_length_slopes(normalized, length)
    association = slope_association(slopes, meta, "condition")
    assert association["p_value"] < 0.01
    log_cpm = np.log2(normalized + 1)
    lfc = (
        log_cpm[meta["condition"] == "case"].mean()
        - log_cpm[meta["condition"] == "control"].mean()
    )
    trend, binned = lfc_length_trend(
        pl.DataFrame({"gene_id": lfc.index, "log2_fold_change": lfc.to_numpy()}), length
    )
    assert trend["spearman_rho"] > 0.5
    assert binned.height == 10


def test_no_length_bias_without_planting():
    counts, meta, length, _ = simulate_bulk(n_de=0, seed=3)
    normalized = counts.div(counts.sum(axis=1) / counts.sum(axis=1).mean(), axis=0)
    association = slope_association(sample_length_slopes(normalized, length), meta, "condition")
    assert association["p_value"] > 0.01


# --- hidden factors
# ---------------------------------------------------------------------------------------------
def test_residual_pcs_recover_hidden_factor_and_k_grid():
    counts, meta, _, hidden = simulate_bulk(hidden=1.0, n_per_group=10, seed=5)
    log_x = np.log2(counts + 1)
    factors = residual_pcs(log_x, meta[["batch", "condition"]], k_max=3)
    assert abs(np.corrcoef(factors["hf1"], hidden)[0, 1]) > 0.9
    run = simple_pipeline(counts)

    def pipeline(m, cols):
        return pl.DataFrame(
            {
                "gene_id": counts.columns,
                "log2_fold_change": np.zeros(counts.shape[1]),
                "p_adjusted": run(m),
            }
        )

    table = k_sensitivity(pipeline, meta, factors, k_grid=[0, 1, 2, 3], k_reference=1)
    assert table["k"].to_list() == [0, 1, 2, 3]
    assert table.filter(pl.col("k") == 1)["jaccard_vs_reference"][0] == 1.0


# --- counts I/O and tximport
# ------------------------------------------------------------------------------------
def test_read_gene_table_strips_versions(tmp_path):
    path = tmp_path / "salmon.merged.gene_counts_length_scaled.tsv"
    path.write_text("gene_id\tgene_name\tS1\tS2\nENSG01.3\tA\t10.4\t3.6\nENSG02.1\tB\t0\t1\n")
    counts, names = read_gene_table(path)
    assert list(counts.columns) == ["ENSG01", "ENSG02"]
    assert names["ENSG01"] == "A"
    assert to_integer_counts(counts).loc["S1", "ENSG01"] == 10
    assert mean_gene_length(counts)["ENSG01"] == pytest.approx(7.0)


@pytest.mark.skipif(shutil.which(RSCRIPT) is None, reason="Rscript with tximport not available")
def test_tximport_lengthscaled_preserves_library_size(tmp_path):
    rows = []
    for sample, reads in (("S1", [100, 50, 30]), ("S2", [10, 80, 60])):
        quant = tmp_path / f"{sample}.sf"
        lengths = [1000, 2000, 1500]
        tpm = np.array(reads) / np.array(lengths)
        tpm = tpm / tpm.sum() * 1e6
        lines = ["Name\tLength\tEffectiveLength\tTPM\tNumReads"]
        lines += [
            f"ENST{i}.1\t{ln}\t{ln - 200}\t{t}\t{r}"
            for i, (ln, t, r) in enumerate(zip(lengths, tpm, reads, strict=True))
        ]
        quant.write_text("\n".join(lines) + "\n")
        rows.append(f"{sample}\t{quant}")
    (tmp_path / "samples.tsv").write_text("sample\tquant\n" + "\n".join(rows) + "\n")
    (tmp_path / "tx2gene.tsv").write_text(
        "transcript_id\tgene_id\nENST0.1\tENSG1.2\nENST1.1\tENSG1.2\nENST2.1\tENSG2.5\n"
    )
    result = subprocess.run(
        [
            RSCRIPT,
            str(TXIMPORT_R),
            "--samples",
            str(tmp_path / "samples.tsv"),
            "--tx2gene",
            str(tmp_path / "tx2gene.tsv"),
            "--out",
            str(tmp_path / "gene"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    counts, _ = read_gene_table(tmp_path / "gene.counts.tsv")
    assert list(counts.columns) == ["ENSG1", "ENSG2"]
    assert counts.sum(axis=1).to_numpy() == pytest.approx([180, 150])


DESEQ2_R = Path(__file__).resolve().parent / "r" / "deseq2_reference.R"


def _r_has(package):
    if shutil.which(RSCRIPT) is None:
        return False
    check = subprocess.run(
        [RSCRIPT, "-e", f"library({package})"], capture_output=True, check=False
    )
    return check.returncode == 0


@pytest.mark.skipif(not _r_has("DESeq2"), reason="Rscript with DESeq2 not available")
def test_pydeseq2_agrees_with_r_deseq2(tmp_path):
    counts, meta, _, _ = simulate_bulk(seed=21)
    genes = filter_genes(counts, meta, "condition")
    counts = counts[genes].astype(int)
    counts.to_csv(tmp_path / "counts.tsv", sep="\t")
    meta.to_csv(tmp_path / "meta.tsv", sep="\t")
    result = subprocess.run(
        [
            RSCRIPT,
            str(DESEQ2_R),
            str(tmp_path / "counts.tsv"),
            str(tmp_path / "meta.tsv"),
            str(tmp_path / "r.tsv"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    r = pl.read_csv(tmp_path / "r.tsv", separator="\t", null_values="NA")
    py, _ = run_deseq2(
        counts, meta, "~batch + condition", ("condition", "case", "control"), shrink=False
    )
    joined = py.join(r, on="gene_id", suffix="_r").drop_nulls(["p_value", "p_value_r"])
    lfc_r = np.corrcoef(joined["log2_fold_change"], joined["log2_fold_change_r"])[0, 1]
    p_rho = stats.spearmanr(joined["p_value"], joined["p_value_r"]).statistic
    found_py = set(py.filter(pl.col("p_adjusted") < 0.05)["gene_id"])
    found_r = set(r.filter(pl.col("p_adjusted") < 0.05)["gene_id"])
    assert lfc_r > 0.99
    assert p_rho > 0.98
    assert len(found_py & found_r) / len(found_py | found_r) > 0.8
