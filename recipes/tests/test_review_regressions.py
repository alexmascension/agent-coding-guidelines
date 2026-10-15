"""Regression tests for the issues found in the external review of 2026-10-14 (each failed before the fix)."""

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pytest
import scanpy as sc
import scipy.sparse as sp

from recipes.bulk.deseq2 import coefficient_name, filter_genes, run_deseq2, validated_counts
from recipes.bulk.fdr_calibration import fdr_calibration
from recipes.bulk.hidden_factors import k_sensitivity, residual_pcs
from recipes.common.gene_ids import annotation_table_from_gtf, resolve_markers
from recipes.common.resources import package_script
from recipes.scrnaseq._clustering import best_match_jaccard
from recipes.scrnaseq.abundance import cell_count_table
from recipes.scrnaseq.doublets import doublet_cluster_review
from recipes.scrnaseq.marker_annotation import annotate_clusters, marker_network
from recipes.scrnaseq.normalization import (
    log_normalize,
    median_size_factors,
    validate_size_factors,
)
from recipes.scrnaseq.pc_robustness import propose_n_pcs
from recipes.scrnaseq.pseudobulk import aggregate_counts, pseudobulk_de
from recipes.scrnaseq.qc import mad_outlier_flags
from recipes.stats.permutation import binomial_upper_bound, permutation_pvalue
from recipes.tests.conftest import simulate_counts


# --- covariate dtypes (P1)
# --------------------------------------------------------------------------------------
def test_numeric_covariate_stays_numeric_in_pseudobulk_design():
    adata = simulate_counts(seed=7)
    ages = {f"s{i}": age for i, age in enumerate([20.0, 30.0, 40.0, 50.0, 60.0, 70.0])}
    adata.obs["age"] = adata.obs["sample"].map(ages).astype(float)
    pseudobulk, _ = aggregate_counts(adata, "sample", "cell_type", covariates=["condition", "age"])
    _, meta = pseudobulk["A"]
    assert pd.api.types.is_float_dtype(meta["age"])
    results, _ = pseudobulk_de(
        pseudobulk, design="~age + condition", contrast=("condition", "case", "control")
    )
    assert results.height > 0


def test_missing_covariate_is_rejected():
    adata = simulate_counts(seed=7)
    adata.obs["rin"] = 7.0
    adata.obs.loc[adata.obs["sample"] == "s0", "rin"] = np.nan
    with pytest.raises(ValueError, match="missing values"):
        aggregate_counts(adata, "sample", "cell_type", covariates=["rin"])


def test_pseudobulk_report_lists_absent_combinations_and_empty_schema():
    adata = simulate_counts(seed=7)
    adata = adata[~((adata.obs["sample"] == "s0") & (adata.obs["cell_type"] == "B"))].copy()
    _, report = aggregate_counts(adata, "sample", "cell_type", covariates=["condition"])
    status = report.filter((pl.col("unit") == "s0") & (pl.col("cell_type") == "B"))["status"][0]
    assert status == "absent"
    results, skipped = pseudobulk_de(
        {}, design="~condition", contrast=("condition", "case", "control")
    )
    assert results.height == 0 and "p_adjusted_global" in results.columns


def test_abundance_manifest_keeps_samples_without_cells():
    adata = simulate_counts(seed=7)
    obs = adata.obs[adata.obs["sample"] != "s5"]
    manifest = adata.obs.drop_duplicates("sample").set_index("sample")[["condition"]]
    table = cell_count_table(obs, "sample", "cell_type", ["condition"], manifest=manifest)
    lost = table.filter(pl.col("sample") == "s5")
    assert set(lost["sample_status"]) == {"no_cells"}
    assert lost["proportion"].is_null().all()


# --- permutation and FDR diagnostic (P1)
# ------------------------------------------------------------------------
def test_nan_statistic_is_rejected():
    with pytest.raises(ValueError, match="not finite"):
        permutation_pvalue(lambda labels: float("nan"), np.array([0, 1] * 5), n_permutations=99)


def test_two_sided_pvalue_uses_scipy_convention():
    rng = np.random.default_rng(0)
    values = (
        rng.normal(size=20) + 5.0
    )  # statistic not centred at 0: absolute values would be wrong
    labels = np.array([0, 1] * 10)
    p, _, _ = permutation_pvalue(lambda lab: values[lab == 1].mean(), labels, n_permutations=499)
    assert p > 0.05


def test_binomial_upper_bound_for_zero_events():
    assert binomial_upper_bound(0, 20) == pytest.approx(1 - 0.05 ** (1 / 20), rel=1e-6)


def test_fdr_diagnostic_rejects_invalid_adjusted_pvalues():
    meta = pd.DataFrame({"condition": ["a", "b"] * 4})
    with pytest.raises(ValueError, match="outside"):
        fdr_calibration(lambda m: np.array([1.5, 0.1]), meta, "condition", n_permutations=2)


def test_fdr_diagnostic_reports_uncertainty():
    meta = pd.DataFrame({"condition": ["a", "b"] * 4})
    _, summary = fdr_calibration(
        lambda m: np.array([0.9, np.nan]), meta, "condition", n_permutations=20
    )
    assert summary["fraction_with_discoveries"] == 0
    assert summary["upper_bound_any_discovery"] == pytest.approx(0.139, abs=0.001)


# --- doublets (P1)
# ----------------------------------------------------------------------------------------------
def _two_lineage_adata(hybrid: bool):
    adata = simulate_counts(seed=11)
    keep = np.flatnonzero(adata.obs["cell_type"].isin(["A", "B"]).to_numpy())
    adata = adata[
        np.r_[keep[:20], keep[-20:]]
    ].copy()  # 20 A + 20 B cells: a small cluster, as doublets are
    if hybrid:  # true doublets: the same cells express markers of A and B
        x = adata.layers["counts"].toarray()
        x[:, 0:20] += 30
        adata.layers["counts"] = sp.csr_matrix(x)
    adata.obs["first_clusters"] = "mixed"
    other = simulate_counts(seed=12)
    other = other[other.obs["cell_type"] == "C"].copy()
    other.obs["first_clusters"] = "C"
    combined = ad.concat([adata, other], index_unique="-")
    combined.obs["doublet_any"] = combined.obs["first_clusters"] == "mixed"
    combined.obs["doublet_consensus"] = False
    return combined


def test_cluster_of_mixed_singlets_is_not_removed():
    adata = _two_lineage_adata(hybrid=False)
    lineages = {"A": list(adata.var_names[0:10]), "B": list(adata.var_names[10:20])}
    review = doublet_cluster_review(adata, "first_clusters", "capture", lineage_markers=lineages)
    row = review.filter(pl.col("cluster") == "mixed").row(0, named=True)
    assert row["frac_coexpressing"] < 0.2
    assert row["status"] != "remove"


def test_cluster_of_true_doublets_is_removed():
    adata = _two_lineage_adata(hybrid=True)
    lineages = {"A": list(adata.var_names[0:10]), "B": list(adata.var_names[10:20])}
    review = doublet_cluster_review(adata, "first_clusters", "capture", lineage_markers=lineages)
    assert review.filter(pl.col("cluster") == "mixed")["status"][0] == "remove"


def test_enrichment_is_adjusted_for_capture_rates():
    n = 400
    obs = pd.DataFrame(
        {
            "capture": ["dense"] * 200 + ["sparse"] * 200,
            "first_clusters": ["k1"] * 100 + ["k2"] * 100 + ["k2"] * 200,
        }
    )
    rng = np.random.default_rng(0)
    calls = np.r_[rng.random(200) < 0.2, rng.random(200) < 0.02]
    obs["doublet_any"], obs["doublet_consensus"] = calls, False
    adata = ad.AnnData(X=sp.csr_matrix((n, 5)), obs=obs.set_axis([f"c{i}" for i in range(n)]))
    adata.layers["counts"] = adata.X
    review = doublet_cluster_review(adata, "first_clusters", "capture")
    k1 = review.filter(pl.col("cluster") == "k1").row(0, named=True)
    assert (
        0.5 < k1["fold_vs_expected"] < 1.5
    )  # k1 sits in the dense capture but is not enriched within it


# --- QC (P1)
# ----------------------------------------------------------------------------------------------------
def _qc_adata(mito, nuclear):
    n = len(mito)
    obs = pd.DataFrame(
        {
            "capture": "c1",
            "log1p_total_counts": np.log1p(np.full(n, 1000.0)) + np.arange(n) * 1e-3,
            "log1p_n_genes": np.log1p(np.full(n, 500.0)) + np.arange(n) * 1e-3,
            "mito_fraction": mito,
            "nuclear_fraction": nuclear,
        },
        index=[f"c{i}" for i in range(n)],
    )
    return ad.AnnData(X=sp.csr_matrix((n, 3)), obs=obs)


def test_mad_zero_and_nan_are_not_silently_passed():
    adata = _qc_adata([0.0] * 9 + [0.9], [0.3] * 9 + [np.nan])
    report = mad_outlier_flags(adata, "capture")
    assert (
        report.filter(pl.col("metric") == "mito_fraction")["status"][0] == "not_evaluated_mad_zero"
    )
    assert adata.obs["qc_invalid"].iloc[-1]
    assert adata.obs["qc_fail"].iloc[-1]


def test_absolute_mito_threshold_applies_when_mad_is_zero():
    adata = _qc_adata([0.0] * 9 + [0.9], [0.3] * 10)
    report = mad_outlier_flags(adata, "capture", max_mito_fraction=0.05)
    assert (
        report.filter(pl.col("metric") == "mito_fraction")["status"][0]
        == "absolute_threshold_only"
    )
    assert adata.obs["qc_outlier_mito_fraction"].iloc[-1]


# --- normalization (P1/P2)
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [np.nan, np.inf, 0.0, -1.0])
def test_invalid_size_factors_are_rejected(bad):
    with pytest.raises(ValueError, match="finite and positive"):
        validate_size_factors(np.array([1.0, bad, 2.0]), 3)


def test_median_normalization_matches_scanpy_normalize_total():
    adata = simulate_counts(seed=3)
    log_normalize(adata, median_size_factors(adata))
    reference = ad.AnnData(X=adata.layers["counts"].copy())
    sc.pp.normalize_total(reference, target_sum=None)
    sc.pp.log1p(reference)
    assert np.allclose(adata.layers["log1p_norm"].toarray(), reference.X.toarray(), atol=1e-5)


def test_r_scripts_are_package_resources():
    assert package_script("recipes.scrnaseq", "scran_size_factors.R").is_file()
    assert package_script("recipes.bulk", "tximport_lengthscaled.R").is_file()


# --- PCs (P2)
# ---------------------------------------------------------------------------------------------------
def test_pc_proposal_is_not_truncated():
    singular = np.r_[np.linspace(100, 50, 60), np.full(140, 1.0)]
    variance = (singular[:50] ** 2) / np.sum(singular**2)
    proposal = propose_n_pcs(variance, singular, shape=(200, 1000))
    assert proposal["gavish_donoho"] == 60
    assert proposal["proposal"] == 60
    assert proposal["status"] == "recompute_pca"


# --- bulk (P2)
# --------------------------------------------------------------------------------------------------
def test_filter_ignores_unused_categories():
    counts = pd.DataFrame({"g1": [0, 0, 0, 0], "g2": [20, 20, 20, 20]}, index=list("wxyz"))
    meta = pd.DataFrame(
        {"condition": pd.Categorical(["a", "a", "b", "b"], categories=["a", "b", "c"])},
        index=list("wxyz"),
    )
    assert list(filter_genes(counts, meta, "condition")) == ["g2"]


def test_coefficient_name_is_exact():
    columns = ["Intercept", "group[T.A]", "group[T.AA]"]
    assert coefficient_name(columns, "group", "A") == "group[T.A]"
    with pytest.raises(ValueError, match="interactions"):
        coefficient_name(["Intercept", "group[T.A]", "group[T.A]:sex[T.m]"], "group", "A")


def test_non_integer_counts_need_explicit_permission():
    counts = pd.DataFrame({"g": [1.5, 2.0]})
    with pytest.raises(ValueError, match="non-integer"):
        validated_counts(counts)
    assert validated_counts(counts, allow_non_integer=True)["g"].tolist() == [2, 2]


def test_run_deseq2_shrunk_lfc_is_finite():
    from recipes.tests.test_bulk_recipes import simulate_bulk

    counts, meta, _, _ = simulate_bulk(seed=8)
    meta["age"] = np.linspace(30, 80, len(meta))
    genes = filter_genes(counts, meta, "condition")
    results, _ = run_deseq2(
        counts[genes], meta, "~age + condition", ("condition", "case", "control")
    )
    assert np.isfinite(results["log2_fold_change_shrunk"].to_numpy()).all()


def test_residual_pcs_respect_rank():
    vst = pd.DataFrame(np.random.default_rng(0).normal(size=(5, 2)), index=list("abcde"))
    design = pd.DataFrame({"condition": list("xxyyy")}, index=list("abcde"))
    with pytest.raises(ValueError, match="exceeds"):
        residual_pcs(vst, design, k_max=3)


def test_k_sensitivity_uses_the_factors():
    rng = np.random.default_rng(1)
    n = 20
    meta = pd.DataFrame({"condition": ["a", "b"] * 10}, index=[f"s{i}" for i in range(n)])
    factor = rng.normal(size=n)
    signal = np.outer(factor, rng.normal(size=50)) * 3 + rng.normal(size=(n, 50))
    factors = pd.DataFrame({"hf1": factor}, index=meta.index)

    def pipeline(m, cols):  # residual variance shrinks only if the factor is actually used
        x = np.column_stack(
            [np.ones(n), (m["condition"] == "b").to_numpy(float), *[m[c] for c in cols]]
        )
        beta, *_ = np.linalg.lstsq(x, signal, rcond=None)
        resid = signal - x @ beta
        p = np.clip(resid.var(axis=0) / 10, 0, 1)
        return pl.DataFrame(
            {"gene_id": [f"g{i}" for i in range(50)], "log2_fold_change": beta[1], "p_adjusted": p}
        )

    table = k_sensitivity(pipeline, meta, factors, k_grid=[0, 1], k_reference=1)
    assert (
        table.filter(pl.col("k") == 0)["n_discoveries"][0]
        < table.filter(pl.col("k") == 1)["n_discoveries"][0]
    )


# --- annotation and markers (P2)
# --------------------------------------------------------------------------------
def test_contradictory_marker_signs_are_rejected():
    with pytest.raises(ValueError, match="both positive and negative"):
        marker_network({"A": ["g1", "g2"]}, negative={"A": ["g2"]})


def test_annotation_handles_ties_unscored_and_alignment():
    scores = pd.DataFrame(
        {"A": [1.0, 1.0, np.nan], "B": [0.5, 1.0, np.nan]}, index=["c1", "c2", "c3"]
    )
    clusters = pd.Series(["k", "k", "u"], index=["c1", "c2", "c3"])
    table = annotate_clusters(scores, clusters)
    assert table.filter(pl.col("cluster") == "u")["label"][0] == "unassigned"
    with pytest.raises(ValueError, match="same cell index"):
        annotate_clusters(scores, clusters.set_axis(["x", "y", "z"]))


def test_unknown_multi_id_rule_is_rejected(gtf_path):
    with pytest.raises(ValueError, match="multi_id_rule"):
        resolve_markers({"x": ["DUP"]}, annotation_table_from_gtf(gtf_path), multi_id_rule="first")


def test_best_match_jaccard_contingency_matches_definition():
    rng = np.random.default_rng(0)
    a, b = rng.integers(0, 5, 300), rng.integers(0, 7, 300)
    expected = {}
    for ka in np.unique(a):
        expected[int(ka)] = max(
            np.sum((a == ka) & (b == kb)) / np.sum((a == ka) | (b == kb)) for kb in np.unique(b)
        )
    assert best_match_jaccard(a, b) == pytest.approx(expected)
