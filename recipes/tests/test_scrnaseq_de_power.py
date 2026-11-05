"""Tests of the recipes for single-cell DE power and pseudoreplication (SCRNASEQ §2, §10; STATISTICS §10)."""

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pytest
import scipy.sparse as sp

from recipes.scrnaseq.depth_matching import match_depth
from recipes.scrnaseq.empirical_power import empirical_power
from recipes.scrnaseq.pseudobulk import (
    aggregate_counts,
    detectability_report,
    min_cells_diagnostic,
    null_split_de,
    power_inputs,
    pseudobulk_de,
)
from recipes.stats.heldout import heldout_splits


def simulate_donors(
    n_per_group=4,
    cells=30,
    n_genes=120,
    donor_sd=0.3,
    de=None,
    cell_types=("A",),
    depth=1.0,
    seed=0,
):
    """Negative-binomial counts for donors of two groups; ``de`` maps gene index -> log2 FC in case.

    ``cells`` and ``depth`` may be dicts by cell type.
    """
    rng = np.random.default_rng(seed)
    base = rng.gamma(2.0, 2.0, size=n_genes) + 0.5
    rows, obs = [], []
    donors = [(f"d{i}", "control") for i in range(n_per_group)] + [
        (f"d{n_per_group + i}", "case") for i in range(n_per_group)
    ]
    for cell_type in cell_types:
        n_cells = cells[cell_type] if isinstance(cells, dict) else cells
        scale = depth[cell_type] if isinstance(depth, dict) else depth
        for donor, condition in donors:
            mean = base * rng.lognormal(0, donor_sd, size=n_genes)
            if de and condition == "case":
                for gene, lfc in de.items():
                    mean[gene] *= 2.0**lfc
            n = n_cells[donor] if isinstance(n_cells, dict) else n_cells
            mu = mean[None, :] * scale * rng.lognormal(0, 0.2, size=(n, 1))
            rows.append(rng.negative_binomial(10, 10 / (10 + mu)))
            obs += [
                {"donor": donor, "sample": donor, "condition": condition, "cell_type": cell_type}
            ] * n
    counts = np.vstack(rows).astype(np.float32)
    obs = pd.DataFrame(obs, index=[f"c{i}" for i in range(counts.shape[0])])
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    adata = ad.AnnData(X=sp.csr_matrix(counts), obs=obs, var=var)
    adata.layers["counts"] = sp.csr_matrix(counts)
    return adata


# --- held-out splits (STAT-RES-008) ---------------------------------------------------------


def test_heldout_splits_are_disjoint_stratified_and_reproducible():
    units = np.array([f"u{i}" for i in range(20)])
    strata = np.array(["a"] * 12 + ["b"] * 8)
    splits = heldout_splits(units, strata, n_splits=5, reference_fraction=0.5, seed=1)
    assert len(splits) == 5
    for reference, validation in splits:
        assert set(reference).isdisjoint(validation)
        assert set(reference) | set(validation) == set(units)
        ref_strata = pd.Series(strata, index=units).loc[reference]
        assert (ref_strata == "a").sum() == 6 and (ref_strata == "b").sum() == 4
    assert len({tuple(sorted(r)) for r, _ in splits}) > 1
    again = heldout_splits(units, strata, n_splits=5, reference_fraction=0.5, seed=1)
    assert all(np.array_equal(a[0], b[0]) for a, b in zip(splits, again, strict=True))


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"units": ["u1", "u1", "u2"]}, "unique"),
        ({"reference_fraction": 1.0}, "reference_fraction"),
        ({"strata": ["a", "a", "b"]}, "fewer than 2"),
    ],
)
def test_heldout_splits_rejects_invalid_input(kwargs, message):
    args = {"units": ["u1", "u2", "u3"], "strata": None, "n_splits": 2}
    args.update(kwargs)
    with pytest.raises(ValueError, match=message):
        heldout_splits(**args)


# --- pseudo-replicates (SC-DE-006) ---------------------------------------------------------


def test_aggregate_counts_rejects_pseudo_replicates(sim_adata):
    rng = np.random.default_rng(0)
    sim_adata.obs["pseudo_unit"] = (
        sim_adata.obs["sample"].astype(str) + "_" + rng.choice(["r1", "r2"], size=sim_adata.n_obs)
    )
    with pytest.raises(ValueError, match="pseudo-replicates"):
        aggregate_counts(
            sim_adata, "pseudo_unit", "cell_type", ["condition"], sample_keys=["sample"]
        )
    pseudobulk, _ = aggregate_counts(
        sim_adata, "sample", "cell_type", ["condition"], sample_keys=["sample"]
    )
    assert set(pseudobulk) == {"A", "B", "C"}


def test_aggregate_counts_accepts_pooled_captures_with_donor_key():
    adata = simulate_donors(n_per_group=2, cells=12, n_genes=20)
    adata.obs["capture"] = np.where(adata.obs["donor"].isin(["d0", "d2"]), "cap1", "cap2")
    pseudobulk, _ = aggregate_counts(
        adata, "donor", "cell_type", ["condition"], sample_keys=["capture", "donor"]
    )
    assert pseudobulk["A"][0].shape[0] == 4
    with pytest.raises(ValueError, match="pseudo-replicates"):
        aggregate_counts(adata, "donor", "cell_type", ["condition"], sample_keys=["capture"])


# --- DE columns and detectability (SC-DE-009, SC-DE-010) ------------------------------------


@pytest.fixture(scope="module")
def de_run():
    from recipes.tests.conftest import simulate_counts

    adata = simulate_counts()
    adata.layers["counts"][:, 299] = 0  # a gene removed by the design-based filter
    adata.layers["counts"].eliminate_zeros()
    pseudobulk, _ = aggregate_counts(adata, "sample", "cell_type", ["condition"])
    counts_a, meta_a = pseudobulk["A"]
    meta_d = meta_a.copy()
    meta_d["condition"] = "control"  # a cell type with only one group: not testable
    pseudobulk["D"] = (counts_a, meta_d)
    contrast = ("condition", "case", "control")
    results, skipped = pseudobulk_de(pseudobulk, "~condition", contrast)
    return adata, pseudobulk, results, skipped, contrast


def test_pseudobulk_de_adds_interval_rank_and_quartile(de_run):
    adata, _, results, skipped, _ = de_run
    assert skipped["cell_type"].to_list() == ["D"]
    a = results.filter(pl.col("cell_type") == "A").to_pandas()
    finite = np.isfinite(a["lfc_se"])
    assert (a.loc[finite, "lfc_ci_low"] < a.loc[finite, "log2_fold_change"]).all()
    assert (a.loc[finite, "lfc_ci_high"] > a.loc[finite, "log2_fold_change"]).all()
    top = a.sort_values("p_value").iloc[0]
    assert top["gene_id"] == adata.var_names[40] and top["rank"] == 1
    significant = a["p_adjusted"] < 0.05
    assert set(a.loc[significant, "significance_quartile"]) <= {1.0, 2.0, 3.0, 4.0}
    assert a.loc[~significant, "significance_quartile"].isna().all()
    assert top["significance_quartile"] == 1


def test_detectability_report_statuses(de_run):
    adata, pseudobulk, results, _, contrast = de_run
    names = adata.var_names
    genes = [names[40], names[100], names[299], "not_a_gene"]
    wide = detectability_report(results, pseudobulk, genes, margin=5.0, contrast=contrast)
    narrow = detectability_report(results, pseudobulk, genes, margin=0.01, contrast=contrast)

    def status(table, cell_type, gene):
        row = table.filter((pl.col("cell_type") == cell_type) & (pl.col("gene_id") == gene))
        return row.row(0, named=True)

    assert status(wide, "A", names[40])["status"] == "changed"
    assert status(wide, "A", names[100])["status"] == "no_change_of_interest"
    assert status(narrow, "A", names[100])["status"] == "inconclusive"
    assert status(wide, "A", names[299])["reason"] == "filtered_low_counts"
    assert status(wide, "A", "not_a_gene")["reason"] == "gene_absent"
    assert status(wide, "D", names[40])["reason"] == "cell_type_not_testable"
    assert status(wide, "A", names[40])["n_units_tested"] == 3
    with pytest.raises(ValueError, match="margin"):
        detectability_report(results, pseudobulk, genes, margin=0.0, contrast=contrast)


# --- global-null diagnostic (SC-DE-007) -----------------------------------------------------


def test_null_split_de_on_null_controls():
    adata = simulate_donors(n_per_group=8, cells=25, n_genes=100, donor_sd=0.3, seed=3)
    pseudobulk, _ = aggregate_counts(adata, "donor", "cell_type", ["condition"])
    contrast = ("condition", "case", "control")
    per_split, summary = null_split_de(pseudobulk, "~condition", contrast, n_splits=4, seed=0)
    row = summary.row(0, named=True)
    assert row["status"] == "tested" and row["n_units"] == 8 and row["n_permutations"] == 4
    assert per_split.height == 4 and (per_split["fraction_label_changed"] > 0).any()
    assert row["median_discoveries"] <= 1
    _, permuted = null_split_de(
        pseudobulk, "~condition", contrast, mode="permute_labels", n_splits=2
    )
    assert permuted.row(0, named=True)["n_units"] == 16
    _, small = null_split_de(
        pseudobulk, "~condition", contrast, n_splits=2, min_samples_per_group=5
    )
    assert small.row(0, named=True)["status"].startswith("not_testable")
    with pytest.raises(ValueError, match="mode"):
        null_split_de(pseudobulk, "~condition", contrast, mode="cells")


# --- minimum cells and power inputs (SC-DE-002, SC-DE-008) ----------------------------------


def test_min_cells_diagnostic_shows_drift_of_low_cell_units():
    cells = {f"d{i}": (3 if i % 4 == 0 else 80) for i in range(24)}
    adata = simulate_donors(n_per_group=12, cells={"A": cells}, n_genes=200, donor_sd=0.1, seed=5)
    pseudobulk, _ = aggregate_counts(adata, "donor", "cell_type", ["condition"], min_cells=1)
    units, by_bin = min_cells_diagnostic(pseudobulk)
    assert units.height == 24
    low = by_bin.filter(pl.col("cells_bin") == "1-4").row(0, named=True)
    high = by_bin.filter(pl.col("cells_bin") == "50-99").row(0, named=True)
    assert low["n_units"] == 6 and high["n_units"] == 18
    assert low["median_distance"] > high["median_distance"]
    assert low["rho_distance_log_cells"] < 0


def test_power_inputs_counts_units_cells_and_umis(sim_adata):
    pseudobulk, _ = aggregate_counts(sim_adata, "sample", "cell_type", ["condition"])
    table = power_inputs(pseudobulk, "condition")
    assert table.height == 6 and set(table["n_units"]) == {3}
    counts, meta = pseudobulk["A"]
    case = meta["condition"].astype(str) == "case"
    expected = np.median(counts[case].sum(axis=1) / meta.loc[case, "n_cells"])
    row = table.filter((pl.col("cell_type") == "A") & (pl.col("group") == "case")).row(
        0, named=True
    )
    assert row["cells_per_unit_median"] == 15
    assert row["umis_per_cell_median"] == pytest.approx(expected)


# --- depth matching (SC-DE-008) -------------------------------------------------------------


def test_match_depth_matches_cells_and_umis():
    adata = simulate_donors(
        n_per_group=3,
        cells={"T": 40, "R": 15},
        n_genes=150,
        cell_types=("T", "R"),
        depth={"T": 2.0, "R": 1.0},
        seed=7,
    )
    matched, per_unit, summary = match_depth(adata, "donor", "cell_type", "T", "R", seed=1)
    assert (per_unit["n_kept"] == 15).all() and matched.n_obs == 90
    assert summary["umi_factor"] < 0.7
    assert summary["mean_umis_target_after"] == pytest.approx(
        summary["mean_umis_reference"], rel=0.03
    )
    original = adata[matched.obs_names].layers["counts"].toarray()
    thinned = matched.layers["counts"].toarray()
    assert np.all(thinned <= original) and np.all(thinned == np.round(thinned))
    again, _, _ = match_depth(adata, "donor", "cell_type", "T", "R", seed=1)
    assert np.array_equal(again.layers["counts"].toarray(), thinned)
    cells_only, _, s_cells = match_depth(adata, "donor", "cell_type", "T", "R", match=("cells",))
    assert s_cells["umi_factor"] == 1.0
    assert np.array_equal(
        cells_only.layers["counts"].toarray(),
        adata[cells_only.obs_names].layers["counts"].toarray(),
    )
    _, _, shallow = match_depth(adata, "donor", "cell_type", "R", "T", match=("umis",))
    assert shallow["umi_factor"] == 1.0
    adata.layers["counts"] = adata.layers["counts"] * 0.5
    with pytest.raises(ValueError, match="integers"):
        match_depth(adata, "donor", "cell_type", "T", "R")


# --- empirical power (SC-DES-007, STAT-RES-008) ---------------------------------------------


def test_empirical_power_recovers_large_effects_and_marks_infeasible_designs():
    de = {**{g: 2.5 for g in range(10)}, **{g: 0.4 for g in range(10, 20)}}
    adata = simulate_donors(n_per_group=12, cells=30, n_genes=120, donor_sd=0.2, de=de, seed=11)
    runs, summary = empirical_power(
        adata,
        "donor",
        "cell_type",
        "A",
        covariates=["condition"],
        design="~condition",
        contrast=("condition", "case", "control"),
        units_per_group=[3, 6, 8],
        cells_per_unit=[10, None],
        n_splits=1,
        n_iterations=1,
        seed=0,
    )
    assert set(runs.filter(pl.col("units_per_group") == 8)["status"]) == {"infeasible"}
    assert set(runs.filter(pl.col("units_per_group") < 8)["status"]) == {"ok"}

    def recall(n_units, n_cells):
        row = summary.filter(
            (pl.col("units_per_group") == n_units)
            & (pl.col("lfc_bin") == "2-3")
            & (
                pl.col("cells_per_unit").is_null()
                if n_cells is None
                else pl.col("cells_per_unit") == n_cells
            )
        )
        return row["mean_recall"][0]

    assert recall(6, None) > 0.5
    assert recall(6, None) >= recall(3, 10)
    assert summary.filter(pl.col("units_per_group") == 8)["n_infeasible"].min() == 1
