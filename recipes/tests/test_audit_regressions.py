"""Regression tests for the external audit of the recipes (C01-C08, C12, C14, C18; layers and PCs).

Each test reproduces the failure the audit reported, or locks a contract stated in the guides
(SCRNASEQ ``SC-DATA-004``, ``SC-PC-002``, ``SC-QC-003``).
"""

import json
import logging
import os
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pytest
import scipy.sparse as sp

import recipes.scrnaseq.pseudobulk as pseudobulk_module
from recipes.bulk.deseq2 import RESULT_COLUMNS, run_deseq2, validated_counts
from recipes.scrnaseq.abundance import cell_count_table
from recipes.scrnaseq.depth_matching import match_depth
from recipes.scrnaseq.doublets import run_scdblfinder
from recipes.scrnaseq.empirical_power import empirical_power
from recipes.scrnaseq.layers import PROTECTED_LAYERS, input_layer, require_integer_counts
from recipes.scrnaseq.normalization import log_normalize, median_size_factors
from recipes.scrnaseq.pc_robustness import elbow_index, propose_n_pcs, robustness_verdict
from recipes.scrnaseq.pseudobulk import aggregate_counts, pseudobulk_de
from recipes.scrnaseq.qc import mad_outlier_flags, qc_metrics
from recipes.tests.conftest import simulate_counts
from recipes.tests.test_scrnaseq_de_power import simulate_donors

REPO_ROOT = Path(__file__).resolve().parents[2]
POWER_KWARGS = {
    "covariates": ["condition"],
    "design": "~condition",
    "contrast": ("condition", "case", "control"),
    "units_per_group": [2],
    "cells_per_unit": [None],
    "n_splits": 1,
    "n_iterations": 1,
}


def _run_python(code: str, **extra_env: str) -> subprocess.CompletedProcess:
    env = {**os.environ, "PYTHONPATH": str(REPO_ROOT), **extra_env}
    return subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, env=env, cwd=REPO_ROOT
    )


def _with_corrected_layer(integer: bool = True) -> ad.AnnData:
    """Simulated counts plus a deliberately different ``counts_ambient_corrected`` layer."""
    adata = simulate_counts(seed=5)
    counts = adata.layers["counts"].toarray()
    corrected = np.floor(counts * 0.5) if integer else counts * 0.9
    adata.layers["counts_ambient_corrected"] = sp.csr_matrix(corrected.astype(np.float32))
    return adata


# --- C01: permutation streams do not depend on PYTHONHASHSEED --------------------------------
def test_permute_within_strata_is_identical_across_hash_seeds():
    code = (
        "import json, numpy as np\n"
        "from recipes.stats.permutation import permute_within_strata\n"
        "labels = np.arange(24)\n"
        "strata = np.array(['batch_x', 'batch_y', 'batch_z', 'batch_w'] * 6)\n"
        "rng = np.random.default_rng(42)\n"
        "out = [permute_within_strata(labels, strata, rng).tolist() for _ in range(3)]\n"
        "print(json.dumps(out))\n"
    )
    outputs = []
    for hash_seed in ("1", "2", "4"):
        result = _run_python(code, PYTHONHASHSEED=hash_seed)
        assert result.returncode == 0, result.stderr
        outputs.append(json.loads(result.stdout))
    assert outputs[0] == outputs[1] == outputs[2]


# --- C02: exact integer check and int64 range -------------------------------------------------
def test_validated_counts_rejects_large_fractional_values():
    with pytest.raises(ValueError, match="non-integer"):
        validated_counts(pd.DataFrame({"g": [100000.4, 3.0]}))
    with pytest.raises(ValueError, match="int64"):
        validated_counts(pd.DataFrame({"g": [2.0**63, 1.0]}))
    large = validated_counts(pd.DataFrame({"g": [2.0**40, 0.0]}))
    assert large["g"].dtype == np.int64 and large["g"].iloc[0] == 2**40


# --- C03: unit identifiers are validated before conversion to strings ------------------------
def test_unit_ids_that_collide_as_strings_are_rejected():
    adata = simulate_counts(seed=7)
    mapping = {"s0": 1, "s1": "1", "s2": "s2", "s3": "s3", "s4": "s4", "s5": "s5"}
    adata.obs["sample"] = adata.obs["sample"].map(mapping).astype(object)
    with pytest.raises(ValueError, match="same identifier"):
        aggregate_counts(adata, "sample", "cell_type", covariates=["condition"])
    with pytest.raises(ValueError, match="same identifier"):
        cell_count_table(adata.obs, "sample", "cell_type", ["condition"])


def test_missing_unit_ids_are_rejected():
    adata = simulate_counts(seed=7)
    adata.obs["sample"] = adata.obs["sample"].astype(object)
    adata.obs.loc[adata.obs.index[0], "sample"] = None
    with pytest.raises(ValueError, match="missing 'sample'"):
        aggregate_counts(adata, "sample", "cell_type", covariates=["condition"])
    with pytest.raises(ValueError, match="missing 'sample'"):
        cell_count_table(adata.obs, "sample", "cell_type", ["condition"])


def test_manifest_index_that_collides_as_strings_is_rejected():
    obs = pd.DataFrame({"sample": ["1"] * 2, "cell_type": ["A", "B"], "condition": ["x"] * 2})
    manifest = pd.DataFrame({"condition": ["x", "x"]}, index=pd.Index([1, "1"], dtype=object))
    with pytest.raises(ValueError, match="same identifier"):
        cell_count_table(obs, "sample", "cell_type", ["condition"], manifest=manifest)


# --- C04: omitted cell types change the denominator only when declared -----------------------
def test_cell_types_that_omit_observed_types_need_an_explicit_subcomposition():
    obs = pd.DataFrame(
        {"sample": ["s1"] * 4, "cell_type": ["A", "B", "B", "B"], "condition": ["case"] * 4}
    )
    with pytest.raises(ValueError, match="omits observed types"):
        cell_count_table(obs, "sample", "cell_type", ["condition"], cell_types=["A"])
    full = cell_count_table(obs, "sample", "cell_type", ["condition"])
    row = full.filter(pl.col("cell_type") == "A").row(0, named=True)
    assert row["proportion"] == pytest.approx(0.25) and row["denominator"] == "all_cell_types"
    sub = cell_count_table(
        obs, "sample", "cell_type", ["condition"], cell_types=["A"], subcomposition=True
    )
    row = sub.row(0, named=True)
    assert (row["n_cells_sample"], row["proportion"]) == (1, 1.0)
    assert row["denominator"] == "listed_cell_types"
    with pytest.raises(ValueError, match="needs the list"):
        cell_count_table(obs, "sample", "cell_type", ["condition"], subcomposition=True)


# --- C05: the manifest does not overwrite obs -------------------------------------------------
def test_manifest_that_contradicts_obs_is_rejected():
    obs = pd.DataFrame(
        {
            "sample": ["s1", "s1", "s2"],
            "cell_type": ["A", "B", "A"],
            "condition": ["case", "case", "control"],
            "age": [30, 30, 40],
        }
    )
    manifest = pd.DataFrame(
        {"condition": ["OTHER", "control", "control"], "age": [30.0, 40.0, 50.0]},
        index=["s1", "s2", "s3"],
    )
    with pytest.raises(ValueError, match="disagree"):
        cell_count_table(obs, "sample", "cell_type", ["condition", "age"], manifest=manifest)
    manifest.loc["s1", "condition"] = "case"  # same values; 30 and 30.0 are the same age
    table = cell_count_table(obs, "sample", "cell_type", ["condition", "age"], manifest=manifest)
    assert set(table["sample"]) == {"s1", "s2", "s3"}


# --- C06 and contract 2: degenerate spectra and criteria evaluated -----------------------------
@pytest.mark.parametrize(
    ("spectrum", "message"),
    [
        ([1.0, 1.0, 1.0, 1.0, 1.0], "constant"),
        ([0.5, np.nan, 0.1], "non-finite"),
        ([0.1, 0.5, 0.2], "order"),
        ([0.5, 0.2], "at least 3"),
    ],
)
def test_elbow_index_rejects_degenerate_spectra(spectrum, message):
    with pytest.raises(ValueError, match=message):
        elbow_index(np.array(spectrum))
    with pytest.raises(ValueError, match=message):
        propose_n_pcs(np.array(spectrum))


def test_propose_n_pcs_records_the_criteria_evaluated(caplog):
    singular = np.r_[np.linspace(100, 50, 10), np.full(40, 1.0)]
    variance = (singular[:30] ** 2) / np.sum(singular**2)
    with caplog.at_level(logging.WARNING, logger="recipes.scrnaseq.pc_robustness"):
        elbow_only = propose_n_pcs(variance)
    assert elbow_only["criteria"] == ["elbow"] and elbow_only["gavish_donoho"] is None
    assert "Gavish-Donoho not evaluated" in caplog.text
    both = propose_n_pcs(variance, singular, shape=(50, 1000))
    assert both["criteria"] == ["elbow", "gavish_donoho"] and both["gavish_donoho"] == 10


# --- C07: clusters of interest must exist ------------------------------------------------------
def test_robustness_verdict_rejects_unknown_clusters_of_interest():
    per_n = pl.DataFrame({"n_pcs": [10, 11], "n_clusters": [1, 1], "ari": [0.95, 0.95]})
    per_cluster = pl.DataFrame(
        {"n_pcs": [10, 11], "ref_cluster": [0, 0], "ref_size": [50, 50], "jaccard": [1.0, 0.1]}
    )
    with pytest.raises(ValueError, match="not reference clusters"):
        robustness_verdict(per_n, per_cluster, clusters_of_interest=[99])
    with pytest.raises(ValueError, match="empty"):
        robustness_verdict(per_n, per_cluster, clusters_of_interest=[])
    verdict = robustness_verdict(per_n, per_cluster, clusters_of_interest=[0])
    assert not verdict["robust"] and verdict["unstable_clusters"] == [0]


# --- C08: derived layers do not overwrite observed counts ---------------------------------------
@pytest.mark.parametrize("key", PROTECTED_LAYERS)
def test_log_normalize_refuses_protected_layers(key):
    adata = simulate_counts(seed=2)
    before = adata.layers["counts"].toarray().copy()
    with pytest.raises(ValueError, match="observed counts"):
        log_normalize(adata, median_size_factors(adata), key_added=key)
    assert np.array_equal(adata.layers["counts"].toarray(), before)


def test_log_normalize_refuses_its_input_layer():
    adata = simulate_counts(seed=2)
    adata.layers["my_counts"] = adata.layers["counts"].copy()
    with pytest.raises(ValueError, match="input layer"):
        log_normalize(adata, median_size_factors(adata), layer="my_counts", key_added="my_counts")


# --- C12: abundance does not need pydeseq2 -----------------------------------------------------
def test_abundance_imports_without_pydeseq2():
    code = (
        "import sys\n"
        "class Block:\n"
        "    def find_spec(self, name, path=None, target=None):\n"
        "        if name == 'pydeseq2' or name.startswith('pydeseq2.'):\n"
        "            raise ModuleNotFoundError(name)\n"
        "        return None\n"
        "sys.meta_path.insert(0, Block())\n"
        "import recipes.scrnaseq.abundance\n"
        "import recipes.scrnaseq.layers\n"
        "assert 'pydeseq2' not in sys.modules\n"
        "print('ok')\n"
    )
    result = _run_python(code)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "ok"


# --- C14: the direct scDblFinder wrapper checks the layer ---------------------------------------
def test_run_scdblfinder_refuses_soupx_layer_before_calling_r():
    adata = _with_corrected_layer()
    adata.uns["ambient_correction"] = {"method": "SoupX", "version": "1.6.2"}
    with pytest.raises(ValueError, match="soupx"):
        run_scdblfinder(
            adata, "capture", layer="counts_ambient_corrected", rscript="no-such-rscript"
        )


# --- contract 1 (SC-DATA-004): input layer per step --------------------------------------------
def test_input_layer_follows_the_table(caplog):
    plain = simulate_counts(seed=1)
    corrected = _with_corrected_layer()
    for step in ("qc", "doublets", "normalization", "pseudobulk"):
        assert input_layer(plain, step) == "counts"
    with caplog.at_level(logging.INFO, logger="recipes.scrnaseq.layers"):
        assert input_layer(corrected, "qc") == "counts"
        assert input_layer(corrected, "doublets") == "counts"
        assert input_layer(corrected, "normalization") == "counts_ambient_corrected"
        assert input_layer(corrected, "pseudobulk") == "counts_ambient_corrected"
    assert "input layer 'counts_ambient_corrected'" in caplog.text
    with pytest.raises(ValueError, match="unknown step"):
        input_layer(plain, "velocity")
    del plain.layers["counts"]
    with pytest.raises(ValueError, match="not in adata.layers"):
        input_layer(plain, "qc")


def test_defaults_read_the_corrected_layer_and_explicit_layers_are_honoured():
    adata = _with_corrected_layer()
    corrected = adata.layers["counts_ambient_corrected"]
    totals = np.asarray(corrected.sum(axis=1)).ravel()
    assert median_size_factors(adata) == pytest.approx(totals / np.median(totals))
    pseudobulk, _ = aggregate_counts(adata, "sample", "cell_type", ["condition"])
    pseudobulk_raw, _ = aggregate_counts(
        adata, "sample", "cell_type", ["condition"], layer="counts"
    )
    total_default = sum(frame.to_numpy().sum() for frame, _ in pseudobulk.values())
    total_raw = sum(frame.to_numpy().sum() for frame, _ in pseudobulk_raw.values())
    assert total_default == pytest.approx(float(corrected.sum()))
    assert total_raw == pytest.approx(float(adata.layers["counts"].sum()))
    qc_metrics(adata)  # QC keeps "counts"
    raw_totals = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
    assert adata.obs["log1p_total_counts"].to_numpy() == pytest.approx(np.log1p(raw_totals))


def test_non_integer_corrected_layer_is_rejected_with_remedies():
    adata = _with_corrected_layer(integer=False)
    for call in (
        lambda: aggregate_counts(adata, "sample", "cell_type", ["condition"]),
        lambda: match_depth(adata, "sample", "cell_type", "A", "B"),
        lambda: empirical_power(adata, "sample", "cell_type", "A", **POWER_KWARGS),
    ):
        with pytest.raises(ValueError) as error:
            call()
        message = str(error.value)
        assert "'counts_ambient_corrected'" in message and "roundToInt" in message
        assert "layer='counts'" in message
    pseudobulk, _ = aggregate_counts(adata, "sample", "cell_type", ["condition"], layer="counts")
    assert pseudobulk
    require_integer_counts(adata, "counts")


# --- G-SC05 (SC-QC-003): fractions are validated before the logit clip --------------------------
def test_mito_percentages_are_refused_not_clipped():
    n = 20
    obs = pd.DataFrame(
        {
            "capture": "c1",
            "log1p_total_counts": np.log1p(np.full(n, 1000.0)) + np.arange(n) * 1e-3,
            "log1p_n_genes": np.log1p(np.full(n, 500.0)) + np.arange(n) * 1e-3,
            "mito_fraction": np.linspace(1.0, 8.0, n),  # percentages
        },
        index=[f"c{i}" for i in range(n)],
    )
    adata = ad.AnnData(X=sp.csr_matrix((n, 3)), obs=obs)
    with pytest.raises(ValueError, match=r"outside \[0, 1\]"):
        mad_outlier_flags(adata, "capture", n_mads_nuclear=None)


# --- C18: empty and non-testable paths have a clear status or error ----------------------------
def test_zero_testable_genes_are_reported_not_testable():
    units = [f"u{i}" for i in range(6)]
    meta = pd.DataFrame({"condition": ["a"] * 3 + ["b"] * 3, "n_cells": [20] * 6}, index=units)
    low = pd.DataFrame(
        np.random.default_rng(0).integers(0, 3, size=(6, 30)),
        index=units,
        columns=[f"g{i}" for i in range(30)],
    )
    results, not_testable = pseudobulk_de(
        {"T": (low, meta)}, "~condition", ("condition", "b", "a")
    )
    assert results.height == 0
    assert not_testable.to_dicts() == [
        {"cell_type": "T", "reason": "no genes pass the gene filter"}
    ]
    with pytest.raises(ValueError, match="no genes"):
        run_deseq2(low.iloc[:, :0], meta, "~condition", ("condition", "b", "a"))


def test_all_non_finite_p_values_leave_global_adjustment_missing(monkeypatch):
    def no_p_values(counts, metadata, design, contrast, **kwargs):
        n = counts.shape[1]
        nan = [np.nan] * n
        table = pl.DataFrame(
            {
                "gene_id": list(counts.columns),
                "base_mean": [1.0] * n,
                "log2_fold_change": [0.0] * n,
                "lfc_se": [1.0] * n,
                "log2_fold_change_shrunk": nan,
                "stat": nan,
                "p_value": nan,
                "p_adjusted": nan,
                "contrast": ["condition_b_vs_a"] * n,
            }
        )
        return table.select(RESULT_COLUMNS), None

    monkeypatch.setattr(pseudobulk_module, "run_deseq2", no_p_values)
    units = [f"u{i}" for i in range(6)]
    meta = pd.DataFrame({"condition": ["a"] * 3 + ["b"] * 3, "n_cells": [20] * 6}, index=units)
    counts = pd.DataFrame(np.full((6, 5), 30), index=units, columns=[f"g{i}" for i in range(5)])
    results, _ = pseudobulk_de({"T": (counts, meta)}, "~condition", ("condition", "b", "a"))
    assert results.height == 5
    assert not np.isfinite(results["p_adjusted_global"].to_numpy()).any()


def test_empirical_power_marks_splits_without_a_testable_reference():
    adata = simulate_donors(n_per_group=6, cells=30, n_genes=60, seed=3)
    runs, summary = empirical_power(
        adata, "donor", "cell_type", "A", min_cells=1000, **POWER_KWARGS
    )
    assert set(runs["status"]) == {"reference_not_testable"}
    assert (summary["n_runs"] == 0).all() and (summary["n_reference_not_testable"] == 1).all()
    small = simulate_donors(n_per_group=4, cells=30, n_genes=60, seed=3)
    runs, _ = empirical_power(
        small, "donor", "cell_type", "A", reference_fraction=0.2, **POWER_KWARGS
    )
    assert set(runs["status"]) == {"reference_not_testable"}
