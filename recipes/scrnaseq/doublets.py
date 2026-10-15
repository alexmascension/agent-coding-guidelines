"""Recipe: doublets with two detectors per capture, consensus/union calls, cluster review.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §4.3 (doublets).

Procedure
---------
1. ``detect_doublets``: per capture, on raw counts, run Scrublet (``scanpy.pp.scrublet``) and
   DoubletDetection
   (``BoostClassifier``, Leiden clustering with standard scaling so that a seed makes it
   reproducible). Store scores,
   calls and the DoubletDetection label (``doublet``, ``singlet``, ``ambiguous``); then
   ``doublet_consensus`` (both
   call doublet) and ``doublet_any`` (either calls doublet). "Not called" is not "confirmed
   singlet": ambiguous
   DoubletDetection cells are kept as ``ambiguous`` in ``dd_label``. The parameters used are
   returned per capture.
2. After a first clustering, ``doublet_cluster_review`` gives per cluster:
   * **capture-adjusted enrichment**: expected ``doublet_any`` calls = sum over captures of (cells
     of the cluster in
     the capture x call rate of the capture); one-sided Poisson test of observed vs expected, BH
     across clusters.
     Rates differ between captures (loading density), so a cluster concentrated in a dense capture
     is not enriched
     for that reason alone. The p-values are exploratory: clusters were defined on the same cells.
   * **per-cell co-expression**: a cell expresses a lineage when its marker score (mean
     log-normalized expression of
     the lineage markers) reaches ``CELL_LEVEL_FRACTION`` of the lineage's reference level: the
     90th percentile of
     its score in its home cluster (the cluster with the highest mean score; the upper quantile
     stays a lineage
     level even when the home cluster mixes lineages). The cluster's co-expression is the fraction
     of its cells expressing two
     or more lineages. A cluster that merely mixes singlets of two lineages has a low fraction.
   Status: ``remove`` (enriched AND co-expression fraction >= ``COEXPRESSION_MIN_FRACTION``),
   ``review`` (one of
   the two), ``keep``.
3. ``doublet_removal_mask``: cells to remove = ``doublet_consensus`` OR cell in a ``remove``
   cluster. ``review``
   clusters are decided manually and the decision is recorded.

The intersection of the two detectors is a conservative project default; no benchmark validates it
as a universal
combination rule (the guide says so).

Defaults (project-adjustable)
-----------------------------
``ENRICHMENT_FDR = 0.05``, ``ENRICHMENT_MIN_FOLD = 2`` (observed at least twice the
capture-adjusted expectation),
``CELL_LEVEL_FRACTION = 0.5``, ``COEXPRESSION_MIN_FRACTION = 0.5``. Scrublet
``expected_doublet_rate`` and the
DoubletDetection thresholds are parameters of ``detect_doublets`` (library defaults unless the
project states the
loading-based expectation).
"""

import logging

import anndata as ad
import doubletdetection
import numpy as np
import polars as pl
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

ENRICHMENT_FDR = 0.05
ENRICHMENT_MIN_FOLD = 2.0
CELL_LEVEL_FRACTION = 0.5
REFERENCE_QUANTILE = (
    0.9  # reference level of a lineage: upper quantile in its home cluster (robust to mixing)
)
COEXPRESSION_MIN_FRACTION = 0.5
DD_N_ITERS = 10  # DoubletDetection fits; its default is 10 in 4.x
DD_P_THRESH = 1e-7  # DoubletDetection 4.x default
DD_VOTER_THRESH = 0.9  # DoubletDetection 4.x default
SCRUBLET_EXPECTED_RATE = (
    0.05  # scanpy default; set it from the loading and the chemistry's documented multiplet rate
)


def detect_doublets(
    adata: ad.AnnData,
    capture_key: str,
    layer: str = "counts",
    seed: int = 0,
    scrublet_expected_rate: float = SCRUBLET_EXPECTED_RATE,
    dd_p_thresh: float = DD_P_THRESH,
    dd_voter_thresh: float = DD_VOTER_THRESH,
) -> pl.DataFrame:
    """Run Scrublet and DoubletDetection per capture and add results to ``obs`` (in place).

    Adds ``scrublet_score``, ``scrublet_call``, ``dd_score``, ``dd_label``, ``dd_call``,
    ``doublet_consensus``,
    ``doublet_any``.

    Parameters
    ----------
    adata : ad.AnnData
        All captures, before integration; raw counts in ``layer``.
    capture_key : str
        ``obs`` column identifying the capture (doublets form within a capture).
    layer : str
        Layer with raw counts.
    seed : int
        Seed for both methods.
    scrublet_expected_rate : float
        Scrublet ``expected_doublet_rate``.
    dd_p_thresh, dd_voter_thresh : float
        DoubletDetection ``predict`` thresholds.

    Returns
    -------
    pl.DataFrame
        Per capture: ``capture``, ``n_cells``, ``n_scrublet``, ``n_dd``, ``n_dd_ambiguous``,
        ``n_consensus``,
        ``n_any``, ``scrublet_threshold``, ``scrublet_expected_rate``, ``dd_p_thresh``,
        ``dd_voter_thresh``.
    """
    n = adata.n_obs
    scores = {k: np.full(n, np.nan) for k in ("scrublet_score", "dd_score")}
    calls = {k: np.zeros(n, dtype=bool) for k in ("scrublet_call", "dd_call")}
    dd_label = np.full(n, "singlet", dtype=object)
    captures = adata.obs[capture_key].astype(str).to_numpy()
    rows = []
    for capture in np.unique(captures):
        idx = np.flatnonzero(captures == capture)
        counts = adata.layers[layer][idx]
        counts = sp.csr_matrix(counts) if not sp.issparse(counts) else counts.tocsr()
        tmp = ad.AnnData(X=counts.astype(np.float32))
        sc.pp.scrublet(
            tmp, expected_doublet_rate=scrublet_expected_rate, random_state=seed, verbose=False
        )
        scores["scrublet_score"][idx] = tmp.obs["doublet_score"].to_numpy()
        calls["scrublet_call"][idx] = tmp.obs["predicted_doublet"].to_numpy(dtype=bool)
        threshold = tmp.uns.get("scrublet", {}).get("threshold", np.nan)

        clf = doubletdetection.BoostClassifier(
            n_iters=DD_N_ITERS,
            clustering_algorithm="leiden",
            standard_scaling=True,
            random_state=seed,
            verbose=False,
        )
        labels = clf.fit(counts).predict(p_thresh=dd_p_thresh, voter_thresh=dd_voter_thresh)
        scores["dd_score"][idx] = clf.doublet_score()
        ambiguous = np.isnan(labels)
        calls["dd_call"][idx] = ~ambiguous & (np.nan_to_num(labels) == 1)
        dd_label[idx[calls["dd_call"][idx]]] = "doublet"
        dd_label[idx[ambiguous]] = "ambiguous"
        rows.append(
            {
                "capture": capture,
                "n_cells": idx.size,
                "n_scrublet": int(calls["scrublet_call"][idx].sum()),
                "n_dd": int(calls["dd_call"][idx].sum()),
                "n_dd_ambiguous": int(ambiguous.sum()),
                "n_consensus": int((calls["scrublet_call"][idx] & calls["dd_call"][idx]).sum()),
                "n_any": int((calls["scrublet_call"][idx] | calls["dd_call"][idx]).sum()),
                "scrublet_threshold": float(threshold),
                "scrublet_expected_rate": scrublet_expected_rate,
                "dd_p_thresh": dd_p_thresh,
                "dd_voter_thresh": dd_voter_thresh,
            }
        )
    for key, values in {**scores, **calls}.items():
        adata.obs[key] = values
    adata.obs["dd_label"] = dd_label
    adata.obs["doublet_consensus"] = calls["scrublet_call"] & calls["dd_call"]
    adata.obs["doublet_any"] = calls["scrublet_call"] | calls["dd_call"]
    return pl.DataFrame(rows)


def combine_calls(call_a: np.ndarray, call_b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (consensus, union) of two call vectors (1 = doublet); NaN is not called."""
    a = np.nan_to_num(np.asarray(call_a, dtype=float), nan=0.0) == 1
    b = np.nan_to_num(np.asarray(call_b, dtype=float), nan=0.0) == 1
    return a & b, a | b


def lineage_scores(
    adata: ad.AnnData, lineage_markers: dict[str, list[str]], layer: str
) -> dict[str, np.ndarray]:
    """Per cell, mean log1p expression of each lineage's markers after median library scaling."""
    x = adata.layers[layer]
    x = sp.csr_matrix(x) if not sp.issparse(x) else x.tocsr()
    totals = np.asarray(x.sum(axis=1)).ravel().astype(float)
    factors = np.where(totals > 0, totals / np.median(totals[totals > 0]), 1.0)
    result = {}
    for lineage, ids in lineage_markers.items():
        cols = np.flatnonzero(adata.var_names.isin(ids))
        if cols.size == 0:
            raise ValueError(f"no markers of lineage {lineage!r} found in var_names")
        result[lineage] = np.log1p(x[:, cols].toarray() / factors[:, None]).mean(axis=1)
    return result


def cells_expressing_lineages(
    scores: dict[str, np.ndarray],
    clusters: np.ndarray,
    cell_level_fraction: float = CELL_LEVEL_FRACTION,
) -> np.ndarray:
    """Count lineages expressed per cell (score >= fraction x 90th percentile in home cluster)."""
    unique = np.unique(clusters)
    n_lineages = np.zeros(clusters.size, dtype=int)
    for values in scores.values():
        means = {c: float(values[clusters == c].mean()) for c in unique}
        home = max(means, key=means.get)
        threshold = cell_level_fraction * float(
            np.quantile(values[clusters == home], REFERENCE_QUANTILE)
        )
        if threshold <= 0:
            continue  # markers undetected even in the home cluster: lineage cannot be scored
        n_lineages += values >= threshold
    return n_lineages


def doublet_cluster_review(
    adata: ad.AnnData,
    cluster_key: str,
    capture_key: str,
    lineage_markers: dict[str, list[str]] | None = None,
    layer: str = "counts",
    fdr: float = ENRICHMENT_FDR,
    min_fold: float = ENRICHMENT_MIN_FOLD,
    cell_level_fraction: float = CELL_LEVEL_FRACTION,
    min_coexpression_fraction: float = COEXPRESSION_MIN_FRACTION,
) -> pl.DataFrame:
    """Tabulate capture-adjusted doublet enrichment and per-cell lineage co-expression per cluster.

    Parameters
    ----------
    adata : ad.AnnData
        Cells with ``doublet_any`` and ``doublet_consensus`` (``detect_doublets``) and a first
        clustering.
    cluster_key : str
        ``obs`` column with the first clustering.
    capture_key : str
        ``obs`` column with the capture; call rates are taken per capture.
    lineage_markers : dict[str, list[str]] | None
        Lineage -> marker IDs (e.g. neuron, astrocyte, oligodendrocyte, microglia), resolved with
        ``recipes.common.gene_ids.resolve_markers``. None skips co-expression (no cluster gets
        ``remove``).
    layer : str
        Raw counts layer.
    fdr, min_fold, cell_level_fraction, min_coexpression_fraction : float
        Thresholds (see module docstring).

    Returns
    -------
    pl.DataFrame
        Per cluster: ``cluster``, ``n_cells``, ``n_any``, ``expected_any``, ``fold_vs_expected``,
        ``frac_consensus``,
        ``p_value``, ``p_adjusted``, ``enriched``, ``frac_coexpressing``, ``coexpressing``,
        ``status``.
    """
    clusters = adata.obs[cluster_key].astype(str).to_numpy()
    captures = adata.obs[capture_key].astype(str).to_numpy()
    any_call = adata.obs["doublet_any"].to_numpy(dtype=bool)
    consensus = adata.obs["doublet_consensus"].to_numpy(dtype=bool)
    capture_rate = {c: float(any_call[captures == c].mean()) for c in np.unique(captures)}
    expected_per_cell = np.array([capture_rate[c] for c in captures])
    n_lineages = None
    if lineage_markers:
        n_lineages = cells_expressing_lineages(
            lineage_scores(adata, lineage_markers, layer), clusters, cell_level_fraction
        )
    rows = []
    for cluster in np.unique(clusters):
        in_c = clusters == cluster
        observed, expected = int(any_call[in_c].sum()), float(expected_per_cell[in_c].sum())
        p = (
            float(poisson.sf(observed - 1, expected))
            if expected > 0
            else (1.0 if observed == 0 else 0.0)
        )
        rows.append(
            {
                "cluster": cluster,
                "n_cells": int(in_c.sum()),
                "n_any": observed,
                "expected_any": expected,
                "fold_vs_expected": observed / expected if expected > 0 else float("nan"),
                "frac_consensus": float(consensus[in_c].mean()),
                "p_value": p,
                "frac_coexpressing": float((n_lineages[in_c] >= 2).mean())
                if n_lineages is not None
                else None,
            }
        )
    table = pl.DataFrame(rows)
    table = table.with_columns(
        p_adjusted=pl.Series(multipletests(table["p_value"].to_numpy(), method="fdr_bh")[1])
    )
    table = table.with_columns(
        enriched=(pl.col("p_adjusted") < fdr)
        & (pl.col("fold_vs_expected").fill_nan(0.0) >= min_fold),
        coexpressing=pl.col("frac_coexpressing").fill_null(0.0) >= min_coexpression_fraction,
    )
    return table.with_columns(
        status=pl.when(pl.col("enriched") & pl.col("coexpressing"))
        .then(pl.lit("remove"))
        .when(pl.col("enriched") | pl.col("coexpressing"))
        .then(pl.lit("review"))
        .otherwise(pl.lit("keep"))
    )


def doublet_removal_mask(adata: ad.AnnData, cluster_key: str, review: pl.DataFrame) -> np.ndarray:
    """Return a mask of cells to remove: ``doublet_consensus`` or in a ``remove`` cluster.

    ``review`` is the output of ``doublet_cluster_review``; set the status of ``review`` clusters
    manually before
    calling (to ``remove`` or ``keep``) and record the decision.
    """
    if (review["status"] == "review").any():
        logger.warning(
            "clusters with status 'review' remain; they are kept unless set to 'remove'"
        )
    remove_clusters = set(review.filter(pl.col("status") == "remove")["cluster"].to_list())
    in_removed = adata.obs[cluster_key].astype(str).isin(remove_clusters).to_numpy()
    return adata.obs["doublet_consensus"].to_numpy(dtype=bool) | in_removed
