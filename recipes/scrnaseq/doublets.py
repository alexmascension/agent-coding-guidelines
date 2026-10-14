"""Recipe: doublets with two detectors per capture, consensus/union calls, and cluster-level review.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §4.3 (doublets).

Procedure
---------
1. ``detect_doublets``: per capture, on raw counts, run Scrublet (``scanpy.pp.scrublet``) and DoubletDetection
   (``BoostClassifier``, Leiden clustering with standard scaling so that a seed makes it reproducible). Store scores
   and calls, then ``doublet_consensus`` (both call doublet) and ``doublet_any`` (either calls doublet).
   DoubletDetection "ambiguous" cells (NaN label) count as not called.
2. After a first clustering, ``doublet_cluster_review`` tabulates per cluster the fraction of ``doublet_any``,
   a one-sided binomial test against the overall rate (Benjamini-Hochberg across clusters), and, if lineage markers
   are given, the two lineages with the highest relative marker level. The relative level of a lineage in a cluster is
   the cluster mean of its marker score (mean log-normalized expression of its markers) divided by the highest cluster
   mean of that lineage, so 1 means "as high as in the lineage's own cluster". Status per cluster:
   ``remove`` (enriched AND co-expressing two lineages), ``review`` (one of the two), ``keep``.
3. ``doublet_removal_mask``: cells to remove = ``doublet_consensus`` OR cell in a ``remove`` cluster. ``review``
   clusters are decided manually and the decision is recorded.

Defaults (project-adjustable)
-----------------------------
``ENRICHMENT_FDR = 0.05``, ``ENRICHMENT_MIN_FOLD = 2`` (cluster rate at least twice the overall rate),
``COEXPRESSION_MIN_LEVEL = 0.5`` (the second lineage reaches at least half of its maximum level).
"""

import logging

import anndata as ad
import doubletdetection
import numpy as np
import polars as pl
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

ENRICHMENT_FDR = 0.05
ENRICHMENT_MIN_FOLD = 2.0
COEXPRESSION_MIN_LEVEL = 0.5
DD_N_ITERS = 10  # DoubletDetection fits; its default is 10 in 4.x


def detect_doublets(adata: ad.AnnData, capture_key: str, layer: str = "counts", seed: int = 0) -> pl.DataFrame:
    """Run Scrublet and DoubletDetection per capture and add calls to ``obs`` (in place).

    Adds ``scrublet_score``, ``scrublet_call``, ``dd_score``, ``dd_call``, ``doublet_consensus``, ``doublet_any``.

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

    Returns
    -------
    pl.DataFrame
        Per capture: ``capture``, ``n_cells``, ``n_scrublet``, ``n_dd``, ``n_consensus``, ``n_any``.
    """
    n = adata.n_obs
    out = {k: np.full(n, np.nan) for k in ("scrublet_score", "dd_score")}
    calls = {k: np.zeros(n, dtype=bool) for k in ("scrublet_call", "dd_call")}
    captures = adata.obs[capture_key].astype(str).to_numpy()
    rows = []
    for capture in np.unique(captures):
        idx = np.flatnonzero(captures == capture)
        counts = adata.layers[layer][idx]
        counts = sp.csr_matrix(counts) if not sp.issparse(counts) else counts.tocsr()
        tmp = ad.AnnData(X=counts.astype(np.float32))
        sc.pp.scrublet(tmp, random_state=seed, verbose=False)
        out["scrublet_score"][idx] = tmp.obs["doublet_score"].to_numpy()
        calls["scrublet_call"][idx] = tmp.obs["predicted_doublet"].to_numpy(dtype=bool)

        clf = doubletdetection.BoostClassifier(
            n_iters=DD_N_ITERS, clustering_algorithm="leiden", standard_scaling=True, random_state=seed, verbose=False
        )
        labels = clf.fit(counts).predict()
        out["dd_score"][idx] = clf.doublet_score()
        calls["dd_call"][idx] = np.nan_to_num(labels, nan=0.0) == 1
        rows.append(
            {
                "capture": capture,
                "n_cells": idx.size,
                "n_scrublet": int(calls["scrublet_call"][idx].sum()),
                "n_dd": int(calls["dd_call"][idx].sum()),
                "n_consensus": int((calls["scrublet_call"][idx] & calls["dd_call"][idx]).sum()),
                "n_any": int((calls["scrublet_call"][idx] | calls["dd_call"][idx]).sum()),
            }
        )
    for key, values in {**out, **calls}.items():
        adata.obs[key] = values
    adata.obs["doublet_consensus"] = calls["scrublet_call"] & calls["dd_call"]
    adata.obs["doublet_any"] = calls["scrublet_call"] | calls["dd_call"]
    return pl.DataFrame(rows)


def combine_calls(call_a: np.ndarray, call_b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (consensus, union) of two boolean call vectors; NaN counts as not called."""
    a = np.nan_to_num(np.asarray(call_a, dtype=float), nan=0.0) == 1
    b = np.nan_to_num(np.asarray(call_b, dtype=float), nan=0.0) == 1
    return a & b, a | b


def _lineage_scores(adata: ad.AnnData, lineage_markers: dict[str, list[str]], layer: str) -> dict[str, np.ndarray]:
    """Per cell, mean log1p expression of each lineage's markers after library-size scaling to the median."""
    x = adata.layers[layer]
    x = sp.csr_matrix(x) if not sp.issparse(x) else x.tocsr()
    totals = np.asarray(x.sum(axis=1)).ravel().astype(float)
    factors = np.where(totals > 0, totals / np.median(totals[totals > 0]), 1.0)
    scores = {}
    for lineage, ids in lineage_markers.items():
        cols = np.flatnonzero(adata.var_names.isin(ids))
        if cols.size == 0:
            raise ValueError(f"no markers of lineage {lineage!r} found in var_names")
        sub = x[:, cols].toarray() / factors[:, None]
        scores[lineage] = np.log1p(sub).mean(axis=1)
    return scores


def doublet_cluster_review(
    adata: ad.AnnData,
    cluster_key: str,
    lineage_markers: dict[str, list[str]] | None = None,
    layer: str = "counts",
    fdr: float = ENRICHMENT_FDR,
    min_fold: float = ENRICHMENT_MIN_FOLD,
    min_coexpression: float = COEXPRESSION_MIN_LEVEL,
) -> pl.DataFrame:
    """Tabulate doublet enrichment and lineage co-expression per cluster of a first clustering.

    Parameters
    ----------
    adata : ad.AnnData
        Cells with ``doublet_any`` and ``doublet_consensus`` (``detect_doublets``) and a first clustering.
    cluster_key : str
        ``obs`` column with the first clustering.
    lineage_markers : dict[str, list[str]] | None
        Lineage -> marker IDs (e.g. neuron, astrocyte, oligodendrocyte, microglia), resolved with
        ``recipes.common.gene_ids.resolve_markers``. None skips co-expression.
    layer : str
        Raw counts layer.
    fdr, min_fold, min_coexpression : float
        Thresholds (see module docstring).

    Returns
    -------
    pl.DataFrame
        Per cluster: ``cluster``, ``n_cells``, ``frac_any``, ``frac_consensus``, ``fold_vs_overall``, ``p_value``,
        ``p_adjusted``, ``enriched``, ``lineage_1``, ``level_lineage_1``, ``lineage_2``, ``level_lineage_2``,
        ``coexpressing``, ``status``.
    """
    clusters = adata.obs[cluster_key].astype(str).to_numpy()
    any_call = adata.obs["doublet_any"].to_numpy(dtype=bool)
    consensus = adata.obs["doublet_consensus"].to_numpy(dtype=bool)
    overall = any_call.mean()
    scores = _lineage_scores(adata, lineage_markers, layer) if lineage_markers else None
    unique_clusters = np.unique(clusters)
    max_level = {}
    if scores is not None:
        for lineage, values in scores.items():
            max_level[lineage] = max(float(values[clusters == c].mean()) for c in unique_clusters)
    rows = []
    for cluster in unique_clusters:
        in_c = clusters == cluster
        k, size = int(any_call[in_c].sum()), int(in_c.sum())
        p = binomtest(k, size, overall, alternative="greater").pvalue if 0 < overall < 1 else 1.0
        row = {
            "cluster": cluster,
            "n_cells": size,
            "frac_any": k / size,
            "frac_consensus": float(consensus[in_c].mean()),
            "fold_vs_overall": (k / size) / overall if overall > 0 else np.nan,
            "p_value": p,
        }
        top = [(None, np.nan), (None, np.nan)]
        if scores is not None:
            levels = (
                (lin, float(v[in_c].mean()) / max_level[lin] if max_level[lin] > 0 else 0.0)
                for lin, v in scores.items()
            )
            ranked = sorted(levels, key=lambda t: -t[1])
            top = (ranked + [(None, np.nan)] * 2)[:2]
        row.update(
            {"lineage_1": top[0][0], "level_lineage_1": top[0][1], "lineage_2": top[1][0], "level_lineage_2": top[1][1]}
        )
        rows.append(row)
    table = pl.DataFrame(rows)
    table = table.with_columns(p_adjusted=pl.Series(multipletests(table["p_value"].to_numpy(), method="fdr_bh")[1]))
    table = table.with_columns(
        enriched=(pl.col("p_adjusted") < fdr) & (pl.col("fold_vs_overall") >= min_fold),
        coexpressing=pl.col("level_lineage_2").fill_nan(0.0).fill_null(0.0) >= min_coexpression,
    )
    table = table.with_columns(
        status=pl.when(pl.col("enriched") & pl.col("coexpressing"))
        .then(pl.lit("remove"))
        .when(pl.col("enriched") | pl.col("coexpressing"))
        .then(pl.lit("review"))
        .otherwise(pl.lit("keep"))
    )
    return table


def doublet_removal_mask(adata: ad.AnnData, cluster_key: str, review: pl.DataFrame) -> np.ndarray:
    """Return a boolean mask of cells to remove: ``doublet_consensus`` or in a cluster with status ``remove``.

    ``review`` is the output of ``doublet_cluster_review``; set the status of ``review`` clusters manually before
    calling (to ``remove`` or ``keep``) and record the decision.
    """
    if (review["status"] == "review").any():
        logger.warning("clusters with status 'review' remain; they are kept unless set to 'remove'")
    remove_clusters = set(review.filter(pl.col("status") == "remove")["cluster"].to_list())
    in_removed = adata.obs[cluster_key].astype(str).isin(remove_clusters).to_numpy()
    return adata.obs["doublet_consensus"].to_numpy(dtype=bool) | in_removed
