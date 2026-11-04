"""Tests for the 4.1.0 rules in recipes/scrnaseq/qc.py: logit mitochondrial rule with a floor,
top-20-gene fraction and flagged fraction per cluster (SC-QC-003, SC-QC-015)."""

import logging

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pytest
import scipy.sparse as sp

from recipes.scrnaseq.qc import (
    flagged_fraction_by_cluster,
    mad_outlier_flags,
    qc_metrics,
)


def _obs_adata(**columns):
    n = len(next(iter(columns.values())))
    obs = pd.DataFrame({"capture": "c1", **columns}, index=[f"c{i}" for i in range(n)])
    return ad.AnnData(X=sp.csr_matrix((n, 3)), obs=obs)


def _mito_adata(extra):
    rng = np.random.default_rng(0)
    base = rng.uniform(0.01, 0.03, 200)
    mito = np.concatenate([base, extra])
    n = mito.size
    return _obs_adata(
        log1p_total_counts=np.log1p(np.full(n, 1000.0)) + rng.normal(0, 0.05, n),
        log1p_n_genes=np.log1p(np.full(n, 500.0)) + rng.normal(0, 0.05, n),
        mito_fraction=mito,
    )


def test_mito_rule_is_upper_tail_on_logit_scale_and_needs_the_floor():
    # 0.0: low mito is not a bad cell; 0.07: beyond the MAD threshold but below the floor; 0.3: damaged
    adata = _mito_adata([0.0, 0.07, 0.3])
    report = mad_outlier_flags(
        adata, "capture", n_mads_nuclear=None, mito_floor=0.08, min_scaled_mad=0.005
    )
    flags = adata.obs["qc_outlier_mito_fraction"].to_numpy()
    assert flags[-3:].tolist() == [False, False, True]
    row = report.filter(pl.col("metric") == "mito_fraction").row(0, named=True)
    assert row["scale"] == "logit"
    assert row["status"] == "evaluated"
    assert row["floor"] == 0.08
    assert 0.03 < row["upper"] < 0.2  # back-transformed to the fraction scale


def test_mito_rule_without_floor_is_reported():
    adata = _mito_adata([0.0, 0.07, 0.3])
    report = mad_outlier_flags(adata, "capture", n_mads_nuclear=None, min_scaled_mad=0.005)
    flags = adata.obs["qc_outlier_mito_fraction"].to_numpy()
    assert flags[-3:].tolist() == [False, True, True]
    status = report.filter(pl.col("metric") == "mito_fraction")["status"][0]
    assert status == "evaluated_no_floor"


def test_top20_fraction_sparse_and_dense_agree():
    rng = np.random.default_rng(1)
    dense = rng.poisson(0.5, size=(30, 60)).astype(np.float32)
    dense[0, :] = 0
    dense[1, :5] = [100, 1, 1, 1, 1]
    expected = np.sort(dense, axis=1)[:, ::-1][:, :20].sum(axis=1) / dense.sum(axis=1).clip(min=1)
    for matrix in (sp.csr_matrix(dense), dense):
        adata = ad.AnnData(X=np.zeros_like(dense))
        adata.layers["counts"] = matrix
        qc_metrics(adata)
        top = adata.obs["top20_fraction"].to_numpy()
        assert np.isnan(top[0])  # no counts: not evaluable, not 0
        np.testing.assert_allclose(top[1:], expected[1:], rtol=1e-6)


def test_top20_fraction_flags_both_tails():
    rng = np.random.default_rng(2)
    top = np.concatenate([rng.normal(0.3, 0.02, 200), [0.95, 0.05]])
    n = top.size
    adata = _obs_adata(
        log1p_total_counts=np.log1p(np.full(n, 1000.0)) + rng.normal(0, 0.05, n),
        log1p_n_genes=np.log1p(np.full(n, 500.0)) + rng.normal(0, 0.05, n),
        top20_fraction=top,
    )
    report = mad_outlier_flags(adata, "capture", n_mads_nuclear=None)
    assert adata.obs["qc_outlier_top20_fraction"].to_numpy()[-2:].all()
    assert report.filter(pl.col("metric") == "top20_fraction")["status"][0] == "evaluated"
    report = mad_outlier_flags(adata, "capture", n_mads_nuclear=None, n_mads_top20=None)
    assert "top20_fraction" not in set(report["metric"])


def test_flagged_fraction_by_cluster_marks_the_cluster_that_loses_most(caplog):
    n = 300
    clusters = np.repeat(["A", "B", "C"], 100)
    flags = np.zeros(n, dtype=bool)
    flags[:50] = True  # A: 50 %
    flags[100:105] = True  # B: 5 %
    flags[200:204] = True  # C: 4 %
    adata = ad.AnnData(
        X=sp.csr_matrix((n, 2)),
        obs=pd.DataFrame(
            {"leiden": clusters, "qc_outlier_mito_fraction": flags, "qc_fail": flags},
            index=[f"c{i}" for i in range(n)],
        ),
    )
    with caplog.at_level(logging.WARNING):
        table = flagged_fraction_by_cluster(adata, "leiden")
    mito = table.filter(pl.col("metric") == "mito_fraction")
    assert mito.filter(pl.col("inspect"))["cluster"].to_list() == ["A"]
    assert mito.filter(pl.col("cluster") == "A")["fraction"][0] == pytest.approx(0.5)
    assert set(table["metric"]) == {"mito_fraction", "qc_fail"}
    assert "inspect before removing" in caplog.text


def test_flagged_fraction_requires_flags():
    adata = ad.AnnData(
        X=sp.csr_matrix((3, 2)),
        obs=pd.DataFrame({"leiden": ["a", "b", "c"]}, index=["x", "y", "z"]),
    )
    with pytest.raises(ValueError, match="mad_outlier_flags"):
        flagged_fraction_by_cluster(adata, "leiden")


def test_minimum_mad_for_fraction_metrics_in_nuclei():
    # nuclei: almost no mitochondrial reads, MAD 0; without a minimum, every nucleus with one read would be flagged
    mito = np.array([0.0] * 190 + [0.01] * 5 + [0.03, 0.2, 0.0, 0.0, 0.0])
    n = mito.size
    rng = np.random.default_rng(3)
    adata = _obs_adata(
        log1p_total_counts=np.log1p(np.full(n, 1000.0)) + rng.normal(0, 0.05, n),
        log1p_n_genes=np.log1p(np.full(n, 500.0)) + rng.normal(0, 0.05, n),
        mito_fraction=mito,
    )
    report = mad_outlier_flags(adata, "capture", n_mads_nuclear=None, n_mads_top20=None)
    row = report.filter(pl.col("metric") == "mito_fraction").row(0, named=True)
    assert row["status"] == "mad_floor_applied"
    assert row["upper"] == pytest.approx(3 * 0.015)
    flagged = adata.obs["qc_outlier_mito_fraction"].to_numpy()
    assert flagged.sum() == 1 and flagged[196]  # only the 20 % nucleus


def test_minimum_mad_on_linear_fraction_and_not_on_log_metrics(caplog):
    n = 50
    adata = _obs_adata(
        log1p_total_counts=np.full(n, np.log1p(1000.0)),  # MAD 0 on a log metric: not evaluated
        log1p_n_genes=np.log1p(np.full(n, 500.0)) + np.linspace(-0.1, 0.1, n),
        top20_fraction=np.full(n, 0.3),
    )
    adata.obs.iloc[-1, adata.obs.columns.get_loc("top20_fraction")] = 0.5
    with caplog.at_level(logging.WARNING):
        report = mad_outlier_flags(adata, "capture", n_mads_nuclear=None)
    status = dict(zip(report["metric"], report["status"], strict=True))
    assert status["log1p_total_counts"] == "not_evaluated_mad_zero"
    assert status["top20_fraction"] == "mad_floor_applied"
    top = adata.obs["qc_outlier_top20_fraction"].to_numpy()
    assert top[-1] and not top[:-1].any()  # 0.5 > 0.3 + 5 x 0.015
    assert "MAD of log1p_total_counts is 0" in caplog.text
