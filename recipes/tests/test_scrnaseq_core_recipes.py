import numpy as np
import pandas as pd
import polars as pl
import pytest

from recipes.common.gene_ids import (
    annotation_table_from_gtf,
    resolve_markers,
    set_var_gene_ids,
    strip_version,
    symbols_for,
)
from recipes.scrnaseq.abundance import cell_count_table, wide_counts
from recipes.scrnaseq.doublets import combine_calls, doublet_cluster_review, doublet_removal_mask
from recipes.scrnaseq.marker_annotation import annotate_clusters, marker_network, score_cells
from recipes.scrnaseq.normalization import log_normalize, median_size_factors
from recipes.scrnaseq.pseudobulk import aggregate_counts, pseudobulk_de
from recipes.scrnaseq.qc import add_nuclear_fraction, flag_gene_groups, mad_outlier_flags, qc_metrics


# --- gene identifiers -------------------------------------------------------------------------------------------
def test_strip_version_keeps_par_suffix():
    assert strip_version("ENSG00000182378.14_PAR_Y") == "ENSG00000182378_PAR_Y"
    assert strip_version("ENSG00000141510") == "ENSG00000141510"


def test_annotation_table_reads_versions_and_missing_symbols(gtf_path):
    table = annotation_table_from_gtf(gtf_path)
    assert table.height == 7
    row = table.filter(pl.col("gene_id") == "ENSG00000000001").row(0, named=True)
    assert row["gene_id_version"] == "ENSG00000000001.3"
    assert table.filter(pl.col("gene_id") == "ENSG00000000007")["gene_symbol"][0] is None


def test_resolve_markers_reports_every_problem(gtf_path):
    table = annotation_table_from_gtf(gtf_path)
    markers = {"astrocyte": ["GFAP", "Aqp4", "DUP", "NOTAGENE"]}
    resolved, report = resolve_markers(markers, table)
    status = dict(report.select("symbol", "status").iter_rows())
    assert status == {"GFAP": "resolved", "Aqp4": "case_mismatch", "DUP": "multiple_ids", "NOTAGENE": "not_found"}
    assert resolved["astrocyte"] == ["ENSG00000000001", "ENSG00000000005", "ENSG00000000006"]


def test_resolve_markers_most_expressed_rule(gtf_path):
    table = annotation_table_from_gtf(gtf_path)
    expression = {"ENSG00000000005": 1.0, "ENSG00000000006": 5.0}
    resolved, _ = resolve_markers({"x": ["DUP"]}, table, multi_id_rule="most_expressed", expression=expression)
    assert resolved["x"] == ["ENSG00000000006"]


def test_symbols_for_falls_back_to_id(gtf_path):
    table = annotation_table_from_gtf(gtf_path)
    assert symbols_for(["ENSG00000000001", "ENSG00000000007"], table) == ["GFAP", "ENSG00000000007"]


# --- QC ---------------------------------------------------------------------------------------------------------
def _annotation_for(adata):
    ids = [strip_version(g) for g in adata.var_names]
    chromosome = ["MT" if 30 <= i < 35 else "1" for i in range(len(ids))]
    return pl.DataFrame(
        {
            "gene_id": ids,
            "gene_id_version": list(adata.var_names),
            "gene_symbol": [f"G{i}" for i in range(len(ids))],
            "chromosome": chromosome,
            "gene_biotype": ["protein_coding"] * len(ids),
        }
    )


def test_mito_flag_comes_from_chromosome_not_prefix(sim_adata):
    annotation = _annotation_for(sim_adata)
    set_var_gene_ids(sim_adata, annotation)
    counts = flag_gene_groups(sim_adata, annotation)
    assert counts["mt"] == 5
    assert sim_adata.var["mt"].to_numpy()[30:35].all()


def test_mad_flags_catch_low_quality_cells(sim_adata):
    annotation = _annotation_for(sim_adata)
    set_var_gene_ids(sim_adata, annotation)
    flag_gene_groups(sim_adata, annotation)
    x = sim_adata.layers["counts"].toarray()
    x[0, :] = 0
    x[0, 30:35] = 500  # a damaged cell: almost only mitochondrial reads
    sim_adata.layers["counts"] = x
    qc_metrics(sim_adata)
    add_nuclear_fraction(sim_adata)
    report = mad_outlier_flags(sim_adata, capture_key="capture")
    assert sim_adata.obs["qc_outlier"].iloc[0]
    assert sim_adata.obs["qc_outlier"].mean() < 0.1
    assert set(report["metric"].unique()) == {
        "log1p_total_counts",
        "log1p_n_genes",
        "mito_fraction",
        "nuclear_fraction",
    }


def test_nuclear_fraction_value(sim_adata):
    add_nuclear_fraction(sim_adata)
    assert np.nanmedian(sim_adata.obs["nuclear_fraction"]) == pytest.approx(0.3, abs=0.05)


# --- doublets ---------------------------------------------------------------------------------------------------
def test_combine_calls_treats_nan_as_not_called():
    consensus, union = combine_calls(np.array([1, 1, 0, np.nan]), np.array([1, 0, 0, 1]))
    assert consensus.tolist() == [True, False, False, False]
    assert union.tolist() == [True, True, False, True]


def test_doublet_cluster_review_and_mask(sim_adata):
    n = sim_adata.n_obs
    rng = np.random.default_rng(0)
    clusters = sim_adata.obs["cell_type"].astype(str).to_numpy().copy()
    hybrid = np.flatnonzero(clusters == "C")[:20]
    clusters[hybrid] = "AB"
    x = sim_adata.layers["counts"].toarray()
    x[hybrid, 0:20] += 30  # express markers of A and B
    sim_adata.layers["counts"] = x
    sim_adata.obs["first_clusters"] = clusters
    any_call = rng.random(n) < 0.03
    any_call[hybrid] = rng.random(hybrid.size) < 0.8
    sim_adata.obs["doublet_any"] = any_call
    sim_adata.obs["doublet_consensus"] = any_call & (rng.random(n) < 0.5)
    lineages = {"A": list(sim_adata.var_names[0:10]), "B": list(sim_adata.var_names[10:20])}
    review = doublet_cluster_review(sim_adata, "first_clusters", lineage_markers=lineages)
    status = dict(review.select("cluster", "status").iter_rows())
    assert status["AB"] == "remove"
    assert status["C"] == "keep"
    mask = doublet_removal_mask(sim_adata, "first_clusters", review)
    assert mask[hybrid].all()


# --- normalization ----------------------------------------------------------------------------------------------
def test_median_size_factors_and_log_normalize(sim_adata):
    factors = median_size_factors(sim_adata)
    assert factors.mean() == pytest.approx(1.0)
    log_normalize(sim_adata, factors)
    expected = np.log1p(sim_adata.layers["counts"][0].toarray().ravel() / sim_adata.obs["size_factor"].iloc[0])
    assert np.allclose(sim_adata.layers["log1p_norm"][0].toarray().ravel(), expected)


# --- marker annotation ------------------------------------------------------------------------------------------
def test_marker_annotation_recovers_types(sim_adata):
    log_normalize(sim_adata, median_size_factors(sim_adata))
    names = sim_adata.var_names
    net = marker_network({"A": list(names[0:10]), "B": list(names[10:20]), "C": list(names[20:30])})
    scores = score_cells(sim_adata, net)
    table = annotate_clusters(scores, sim_adata.obs["cell_type"])
    assert dict(table.select("cluster", "label").iter_rows()) == {"A": "A", "B": "B", "C": "C"}


# --- pseudobulk and abundance -----------------------------------------------------------------------------------
def test_pseudobulk_detects_planted_de_gene(sim_adata):
    pseudobulk, report = aggregate_counts(sim_adata, "sample", "cell_type", covariates=["condition"])
    assert report["kept"].all()
    results, skipped = pseudobulk_de(pseudobulk, design="~condition", contrast=("condition", "case", "control"))
    top = results.filter(pl.col("cell_type") == "A").sort("p_value").row(0, named=True)
    assert top["gene_id"] == sim_adata.var_names[40]
    assert top["log2_fold_change"] > 1
    assert skipped.height == 0


def test_pseudobulk_rejects_inconsistent_covariate(sim_adata):
    sim_adata.obs.loc[sim_adata.obs.index[0], "condition"] = (
        "case" if sim_adata.obs["condition"].iloc[0] == "control" else "control"
    )
    with pytest.raises(ValueError, match="not constant"):
        aggregate_counts(sim_adata, "sample", "cell_type", covariates=["condition"])


def test_cell_count_table_has_explicit_zeros(sim_adata):
    obs = sim_adata.obs[~((sim_adata.obs["sample"] == "s0") & (sim_adata.obs["cell_type"] == "B"))]
    table = cell_count_table(obs, "sample", "cell_type", ["condition"])
    assert table.height == 6 * 3
    zero = table.filter((pl.col("sample") == "s0") & (pl.col("cell_type") == "B"))
    assert zero["n_cells"][0] == 0
    assert table.group_by("sample").agg(pl.col("proportion").sum())["proportion"].to_numpy() == pytest.approx(1.0)
    assert isinstance(wide_counts(table), pd.DataFrame)


def test_detect_doublets_finds_planted_doublets():
    import anndata as ad
    import scipy.sparse as sp

    from recipes.scrnaseq.doublets import detect_doublets
    from recipes.tests.conftest import simulate_counts

    adata = simulate_counts(seed=3, cells_per_type=600)
    adata = adata[adata.obs["capture"] == "s0"].copy()
    rng = np.random.default_rng(0)
    x = adata.layers["counts"].toarray()
    types = adata.obs["cell_type"].to_numpy()
    a, b = rng.choice(np.flatnonzero(types == "A"), 15), rng.choice(np.flatnonzero(types == "B"), 15)
    doublets = x[a] + x[b]
    combined = ad.AnnData(X=sp.csr_matrix(np.vstack([x, doublets])))
    combined.layers["counts"] = combined.X.copy()
    combined.obs["capture"] = "s0"
    combined.obs["planted"] = np.r_[np.zeros(x.shape[0], bool), np.ones(15, bool)]
    summary = detect_doublets(combined, capture_key="capture", seed=0)
    assert summary["n_any"][0] >= summary["n_consensus"][0]
    planted_any = combined.obs.loc[combined.obs["planted"], "doublet_any"].mean()
    singlet_any = combined.obs.loc[~combined.obs["planted"], "doublet_any"].mean()
    assert planted_any > singlet_any


RSCRIPT = __import__("os").environ.get("RECIPES_RSCRIPT", "Rscript")


@pytest.mark.skipif(__import__("shutil").which(RSCRIPT) is None, reason="Rscript with scran not available")
def test_scran_size_factors_track_depth(sim_adata):
    from recipes.scrnaseq.normalization import scran_size_factors

    table = scran_size_factors(sim_adata, cluster_key="cell_type", rscript=RSCRIPT)
    factors = table["size_factor"].to_numpy()
    totals = np.asarray(sim_adata.layers["counts"].sum(axis=1)).ravel()
    assert factors.mean() == pytest.approx(1.0)
    assert np.corrcoef(np.log(factors), np.log(totals))[0, 1] > 0.8


def _composition_sim(seed=0, n_genes=2000, n_per=500):
    """Three cell types with different RNA composition and known per-cell depth factors."""
    import anndata as ad
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    base = rng.gamma(0.6, 1.5, size=n_genes)
    blocks, depth, types = [], [], []
    for t in range(3):
        mu = base.copy()
        mu[t * 100 : (t + 1) * 100] *= 20
        d = rng.lognormal(0, 0.5, size=n_per)
        blocks.append(rng.negative_binomial(5, 5 / (5 + mu[None, :] * d[:, None])))
        depth += list(d)
        types += [str(t)] * n_per
    adata = ad.AnnData(
        X=sp.csr_matrix(np.vstack(blocks).astype(np.float64)),
        obs=pd.DataFrame({"cell_type": types}, index=[f"c{i}" for i in range(len(types))]),
    )
    adata.layers["counts"] = adata.X.copy()
    return adata, np.array(depth) / np.mean(depth)


def test_scranpy_size_factors_beat_library_size():
    from recipes.scrnaseq.normalization import scranpy_size_factors

    adata, truth = _composition_sim()
    factors = scranpy_size_factors(adata, cluster_key="cell_type")
    library = median_size_factors(adata)
    err_scran = np.median(np.abs(np.log(factors) - np.log(truth)))
    err_library = np.median(np.abs(np.log(library) - np.log(truth)))
    assert factors.mean() == pytest.approx(1.0)
    assert err_scran < err_library


def test_scranpy_refuses_huge_dense_copy(sim_adata):
    from recipes.scrnaseq.normalization import scranpy_size_factors

    with pytest.raises(MemoryError):
        scranpy_size_factors(sim_adata, max_dense_gb=1e-6)


@pytest.mark.skipif(__import__("shutil").which(RSCRIPT) is None, reason="Rscript with scran not available")
def test_scranpy_agrees_with_r_scran():
    from recipes.scrnaseq.normalization import scran_size_factors, scranpy_size_factors

    adata, _ = _composition_sim(n_per=300)
    py = scranpy_size_factors(adata, cluster_key="cell_type")
    r = scran_size_factors(adata, cluster_key="cell_type", rscript=RSCRIPT)["size_factor"].to_numpy()
    assert np.corrcoef(np.log(py), np.log(r))[0, 1] > 0.98
