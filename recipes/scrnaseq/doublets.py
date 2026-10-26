"""Recipe: doublets with two or three detectors per capture, consensus/union calls, cluster review.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §4.3 (doublets).

Procedure
---------
1. ``detect_doublets``: per capture, on counts without ambient correction (``SC-QC-012``), run the
   detectors named in ``detectors``: Scrublet (``scanpy.pp.scrublet``) and scDblFinder (R, through
   ``scdblfinder.R``; doublets simulated within each capture) by default, and optionally
   DoubletDetection (``BoostClassifier``, Leiden clustering with standard scaling so that a seed
   makes it reproducible). Store scores and calls, ``doublet_n_evaluated`` and ``doublet_n_calls``,
   then ``doublet_consensus`` (every evaluated detector calls doublet, at least two evaluated) and
   ``doublet_any`` (at least one calls doublet). "Not called" is not "confirmed singlet":
   ambiguous DoubletDetection cells are kept as ``ambiguous`` in ``dd_label``. The parameters used
   are returned per capture.
   Scrublet thresholds (``SC-QC-013``): a value in ``scrublet_thresholds`` is used for its capture;
   otherwise scanpy sets it at the minimum between the modes of the simulated-score histogram.
   When that fails, scanpy calls no doublet; the recipe raises (``on_scrublet_failure="raise"``)
   or marks Scrublet as not evaluated for that capture (``"not_evaluated"``), never as zero
   doublets. Simulated scores and thresholds are kept in ``uns["doublets"]["scrublet"]`` for the
   histogram check.
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
   cluster OR, in pooled captures, a genotype doublet (``genotype_key``, ``SC-QC-014``).
   ``review`` clusters are decided manually and the decision is recorded.

Input layer
-----------
``ambient_method_of_layer`` checks the layer before ``detect_doublets`` and the co-expression
part of ``doublet_cluster_review`` read it. ``counts`` is accepted. ``counts_ambient_corrected`` is
refused unless ``uns["ambient_correction"]["method"]`` is CellBender: SoupX and decontX corrections
leave more doublets undetected, CellBender was the exception [Meier 2026]. Other layer names are
accepted with a warning, since their correction cannot be checked. ``doublet_cluster_review`` also
warns when a SoupX or decontX layer already exists: with these methods the first clustering is
computed on ``counts`` (in the preprint, doublet clusters merged into singlet clusters after
decontX).

The intersection of the detectors is a conservative project default; the benchmark of
[Neavin 2024] recommends calling a singlet only when more than half of the methods do (for two
methods, the union), and the guide says so.

Defaults (project-adjustable)
-----------------------------
``ENRICHMENT_FDR = 0.05``, ``ENRICHMENT_MIN_FOLD = 2`` (observed at least twice the
capture-adjusted expectation),
``CELL_LEVEL_FRACTION = 0.5``, ``COEXPRESSION_MIN_FRACTION = 0.5``. Scrublet
``expected_doublet_rate``, the scDblFinder ``dbr`` and the DoubletDetection thresholds are
parameters of ``detect_doublets`` (library defaults unless the project states the loading-based
expectation).
"""

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import polars as pl
import scanpy as sc
import scipy.io
import scipy.sparse as sp
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests

from recipes.common.resources import package_script

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
DETECTORS = ("scrublet", "scdblfinder", "doubletdetection")
DEFAULT_DETECTORS = ("scrublet", "scdblfinder")  # SC-QC-009; DoubletDetection is an optional third
OBS_PREFIX = {"scrublet": "scrublet", "scdblfinder": "scdblfinder", "doubletdetection": "dd"}
SCRUBLET_FAILURE_POLICIES = ("raise", "not_evaluated")
UNCORRECTED_LAYER = "counts"
AMBIENT_LAYER = "counts_ambient_corrected"
# Ambient methods whose corrected counts hide doublets, and the exception (SC-QC-012, Meier 2026).
MASKING_AMBIENT_METHODS = frozenset({"soupx", "decontx"})
ACCEPTED_AMBIENT_METHODS = frozenset({"cellbender"})


def ambient_method_of_layer(adata: ad.AnnData, layer: str) -> str:
    """Return the ambient correction behind ``layer``; refuse corrections that hide doublets.

    Parameters
    ----------
    adata : ad.AnnData
        Data with ``layer``; ``uns["ambient_correction"]["method"]`` names the method that produced
        ``counts_ambient_corrected`` (``SC-DATA-001``).
    layer : str
        Layer that doublet detection or the co-expression review would read.

    Returns
    -------
    str
        ``"none"`` for ``counts``; the recorded method, lower case, for
        ``counts_ambient_corrected``; ``"unknown"`` for any other layer (with a warning).

    Raises
    ------
    ValueError
        ``layer`` is absent; or it is ``counts_ambient_corrected`` and the recorded method is
        SoupX, decontX or not an accepted method (including no record).
    """
    if layer not in adata.layers:
        raise ValueError(f"layer {layer!r} not in adata.layers")
    if layer == UNCORRECTED_LAYER:
        return "none"
    if layer != AMBIENT_LAYER:
        logger.warning(
            f"layer {layer!r}: its ambient correction cannot be checked; doublets are detected on "
            f"{UNCORRECTED_LAYER!r} (SC-QC-012)"
        )
        return "unknown"
    record = adata.uns.get("ambient_correction") or {}
    method = str(record.get("method", "")).strip().lower()
    if method in MASKING_AMBIENT_METHODS:
        raise ValueError(
            f"{AMBIENT_LAYER!r} was corrected with {method}, which leaves doublets undetected; "
            f"use layer={UNCORRECTED_LAYER!r} (SC-QC-012)"
        )
    if method not in ACCEPTED_AMBIENT_METHODS:
        recorded = method or "no method"
        raise ValueError(
            f"uns['ambient_correction'] records {recorded!r} for {AMBIENT_LAYER!r}; use "
            f"layer={UNCORRECTED_LAYER!r}, or record CellBender if it produced the layer "
            "(SC-DATA-001, SC-QC-012)"
        )
    logger.warning(
        f"doublets on {AMBIENT_LAYER!r} ({method}): a deviation from SC-QC-012; record it in "
        "PROJECT_CONTEXT.md"
    )
    return method


def _warn_if_masking_layer_exists(adata: ad.AnnData) -> None:
    """Warn when a SoupX or decontX layer exists before the doublet review (``SC-QC-012``)."""
    method = str((adata.uns.get("ambient_correction") or {}).get("method", "")).strip().lower()
    if AMBIENT_LAYER in adata.layers and method in MASKING_AMBIENT_METHODS:
        logger.warning(
            f"{AMBIENT_LAYER!r} ({method}) exists before the doublet review: compute the first "
            "clustering on 'counts', since doublet clusters can merge into singlet clusters after "
            "this correction (SC-QC-012)"
        )


def _run_scrublet(
    counts: sp.csr_matrix, expected_rate: float, threshold: float | None, seed: int
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Run Scrublet on one capture; return scores, calls and threshold information.

    ``info["threshold_source"]`` is ``"manual"`` (``threshold`` given), ``"automatic"`` or
    ``"failed"`` (scanpy found no minimum between two modes and called no doublet).
    """
    tmp = ad.AnnData(X=counts.astype(np.float32))
    sc.pp.scrublet(
        tmp,
        expected_doublet_rate=expected_rate,
        threshold=threshold,
        random_state=seed,
        verbose=False,
    )
    uns = tmp.uns.get("scrublet", {})
    used = threshold if threshold is not None else uns.get("threshold")
    if threshold is not None:
        source = "manual"
    elif used is None:
        source = "failed"
    else:
        source = "automatic"
    calls = (
        tmp.obs["predicted_doublet"].to_numpy(dtype=bool)
        if source != "failed"
        else np.zeros(tmp.n_obs, dtype=bool)
    )
    info = {
        "doublet_scores_sim": np.asarray(uns.get("doublet_scores_sim", []), dtype=float),
        "threshold": float(used) if used is not None else float("nan"),
        "threshold_source": source,
    }
    return tmp.obs["doublet_score"].to_numpy(dtype=float), calls, info


def _run_doubletdetection(
    counts: sp.csr_matrix, seed: int, p_thresh: float, voter_thresh: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run DoubletDetection on one capture; return scores, calls and the ambiguous mask."""
    import doubletdetection  # optional third detector (SC-QC-009)

    clf = doubletdetection.BoostClassifier(
        n_iters=DD_N_ITERS,
        clustering_algorithm="leiden",
        standard_scaling=True,
        random_state=seed,
        verbose=False,
    )
    labels = clf.fit(counts).predict(p_thresh=p_thresh, voter_thresh=voter_thresh)
    ambiguous = np.isnan(labels)
    calls = ~ambiguous & (np.nan_to_num(labels) == 1)
    return np.asarray(clf.doublet_score(), dtype=float), calls, ambiguous


def run_scdblfinder(
    adata: ad.AnnData,
    capture_key: str,
    layer: str = "counts",
    seed: int = 0,
    dbr: float | None = None,
    rscript: str = "Rscript",
    workdir: str | Path | None = None,
) -> pl.DataFrame:
    """Doublet scores and calls from scDblFinder through ``scdblfinder.R``.

    Doublets are simulated within each capture (``samples`` of ``scDblFinder``).

    Parameters
    ----------
    adata : ad.AnnData
        Cells of all captures; counts without ambient correction in ``layer``.
    capture_key : str
        ``obs`` column identifying the capture.
    layer : str
        Layer with the counts.
    seed : int
        Seed set before ``scDblFinder``.
    dbr : float | None
        Expected doublet rate; None keeps the scDblFinder default.
    rscript : str
        ``Rscript`` of an environment with scDblFinder, SingleCellExperiment, Matrix and optparse.
    workdir : str | Path | None
        Parent directory for the exchange files (a unique subdirectory is created and kept, with
        the R provenance file); None uses a temporary directory deleted afterwards.

    Returns
    -------
    pl.DataFrame
        ``barcode``, ``capture``, ``score``, ``class`` (``singlet`` or ``doublet``), in the order
        of ``adata.obs_names``.

    Raises
    ------
    RuntimeError
        If ``rscript`` is missing, the R script fails (its stderr is included) or the barcode order
        returned differs.
    """
    if shutil.which(rscript) is None:
        raise RuntimeError(f"{rscript} not found; run in an environment with R and scDblFinder")
    if workdir is not None:
        Path(workdir).mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="scdblfinder_", dir=workdir))
    r_script = package_script("recipes.scrnaseq", "scdblfinder.R")
    try:
        counts = adata.layers[layer]
        counts = counts if sp.issparse(counts) else sp.csr_matrix(counts)
        scipy.io.mmwrite(tmp / "counts.mtx", counts.T.tocoo())
        (tmp / "barcodes.txt").write_text("\n".join(adata.obs_names) + "\n")
        (tmp / "genes.txt").write_text("\n".join(adata.var_names) + "\n")
        (tmp / "samples.txt").write_text("\n".join(adata.obs[capture_key].astype(str)) + "\n")
        command = [
            rscript,
            str(r_script),
            "--counts",
            str(tmp / "counts.mtx"),
            "--barcodes",
            str(tmp / "barcodes.txt"),
            "--genes",
            str(tmp / "genes.txt"),
            "--samples",
            str(tmp / "samples.txt"),
            "--seed",
            str(seed),
            "--out",
            str(tmp / "scdblfinder.tsv"),
        ]
        if dbr is not None:
            command += ["--dbr", str(dbr)]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"scdblfinder.R failed:\n{result.stderr}")
        table = pl.read_csv(
            tmp / "scdblfinder.tsv",
            separator="\t",
            schema_overrides={"barcode": pl.Utf8, "capture": pl.Utf8, "class": pl.Utf8},
        )
    finally:
        if workdir is None:
            shutil.rmtree(tmp, ignore_errors=True)
    if table["barcode"].to_list() != list(adata.obs_names):
        raise RuntimeError("barcode order returned by R differs from adata.obs_names")
    return table


def detect_doublets(
    adata: ad.AnnData,
    capture_key: str,
    layer: str = "counts",
    detectors: tuple[str, ...] = DEFAULT_DETECTORS,
    seed: int = 0,
    scrublet_expected_rate: float = SCRUBLET_EXPECTED_RATE,
    scrublet_thresholds: dict[str, float] | None = None,
    on_scrublet_failure: str = "raise",
    scdblfinder_dbr: float | None = None,
    dd_p_thresh: float = DD_P_THRESH,
    dd_voter_thresh: float = DD_VOTER_THRESH,
    rscript: str = "Rscript",
    workdir: str | Path | None = None,
) -> pl.DataFrame:
    """Run the doublet detectors per capture and add results to ``obs`` (in place).

    Adds, per detector, ``<prefix>_score`` and ``<prefix>_call`` (prefixes ``scrublet``,
    ``scdblfinder``, ``dd``), ``scrublet_evaluated`` and ``dd_label`` when those detectors run,
    and ``doublet_n_evaluated``, ``doublet_n_calls``, ``doublet_consensus``, ``doublet_any``.
    Scrublet simulated scores and thresholds go to ``uns["doublets"]["scrublet"][capture]``.

    Parameters
    ----------
    adata : ad.AnnData
        All captures, before integration; counts in ``layer``.
    capture_key : str
        ``obs`` column identifying the capture (doublets form within a capture).
    layer : str
        Layer with counts without ambient correction (``counts``); checked with
        ``ambient_method_of_layer``.
    detectors : tuple[str, ...]
        At least two of ``DETECTORS``; default Scrublet and scDblFinder (``SC-QC-009``).
    seed : int
        Seed for every detector.
    scrublet_expected_rate : float
        Scrublet ``expected_doublet_rate``.
    scrublet_thresholds : dict[str, float] | None
        Manual Scrublet thresholds by capture (``SC-QC-013``); other captures use the automatic
        threshold.
    on_scrublet_failure : str
        ``"raise"`` (default) or ``"not_evaluated"`` when the automatic threshold is not found.
    scdblfinder_dbr : float | None
        scDblFinder expected doublet rate; None keeps its default.
    dd_p_thresh, dd_voter_thresh : float
        DoubletDetection ``predict`` thresholds.
    rscript : str
        ``Rscript`` used for scDblFinder.
    workdir : str | Path | None
        Parent directory for the scDblFinder exchange files (see ``run_scdblfinder``).

    Returns
    -------
    pl.DataFrame
        Per capture: ``capture``, ``n_cells``, ``input_layer``, ``ambient_method``,
        ``detectors``, ``n_<prefix>`` per detector, ``n_consensus``, ``n_any``, and the parameters
        of the detectors run (``scrublet_threshold``, ``scrublet_threshold_source``,
        ``scrublet_expected_rate``, ``scdblfinder_dbr``, ``n_dd_ambiguous``, ``dd_p_thresh``,
        ``dd_voter_thresh``).

    Raises
    ------
    ValueError
        Fewer than two or unknown detectors; thresholds for unknown captures; ``layer`` holds
        SoupX- or decontX-corrected counts or an unrecorded correction
        (``ambient_method_of_layer``); an automatic Scrublet threshold failed and
        ``on_scrublet_failure="raise"`` (the diagnostics are stored first).
    """
    detectors = tuple(dict.fromkeys(detectors))
    unknown = sorted(set(detectors) - set(DETECTORS))
    if unknown:
        raise ValueError(f"unknown detectors {unknown}; use some of {DETECTORS}")
    if len(detectors) < 2:
        raise ValueError(f"SC-QC-009 asks for at least two detectors, got {detectors}")
    if on_scrublet_failure not in SCRUBLET_FAILURE_POLICIES:
        raise ValueError(f"on_scrublet_failure must be one of {SCRUBLET_FAILURE_POLICIES}")
    ambient_method = ambient_method_of_layer(adata, layer)
    captures = adata.obs[capture_key].astype(str).to_numpy()
    thresholds = {str(k): float(v) for k, v in (scrublet_thresholds or {}).items()}
    unknown_captures = sorted(set(thresholds) - set(captures))
    if unknown_captures:
        raise ValueError(f"scrublet_thresholds for unknown captures: {unknown_captures}")
    n = adata.n_obs
    score = {det: np.full(n, np.nan) for det in detectors}
    call = {det: np.zeros(n, dtype=bool) for det in detectors}
    evaluated = {det: np.ones(n, dtype=bool) for det in detectors}
    dd_label = np.full(n, "singlet", dtype=object)
    scrublet_uns, failed, rows = {}, [], {}
    for capture in np.unique(captures):
        idx = np.flatnonzero(captures == capture)
        counts = adata.layers[layer][idx]
        counts = sp.csr_matrix(counts) if not sp.issparse(counts) else counts.tocsr()
        row = {"capture": capture, "n_cells": idx.size, "input_layer": layer}
        row |= {"ambient_method": ambient_method, "detectors": ",".join(detectors)}
        if "scrublet" in detectors:
            s_score, s_call, info = _run_scrublet(
                counts, scrublet_expected_rate, thresholds.get(capture), seed
            )
            score["scrublet"][idx], call["scrublet"][idx] = s_score, s_call
            if info["threshold_source"] == "failed":
                failed.append(capture)
                evaluated["scrublet"][idx] = False
            scrublet_uns[capture] = info
            row |= {
                "scrublet_threshold": info["threshold"],
                "scrublet_threshold_source": info["threshold_source"],
                "scrublet_expected_rate": scrublet_expected_rate,
            }
        if "doubletdetection" in detectors:
            d_score, d_call, ambiguous = _run_doubletdetection(
                counts, seed, dd_p_thresh, dd_voter_thresh
            )
            score["doubletdetection"][idx], call["doubletdetection"][idx] = d_score, d_call
            dd_label[idx[d_call]] = "doublet"
            dd_label[idx[ambiguous]] = "ambiguous"
            row |= {
                "n_dd_ambiguous": int(ambiguous.sum()),
                "dd_p_thresh": dd_p_thresh,
                "dd_voter_thresh": dd_voter_thresh,
            }
        rows[capture] = row
    if "scrublet" in detectors:
        adata.uns.setdefault("doublets", {})["scrublet"] = scrublet_uns
    if failed:
        message = (
            f"automatic Scrublet threshold not found in captures {failed}: inspect the histograms "
            "in uns['doublets']['scrublet'] and pass scrublet_thresholds, or use "
            "on_scrublet_failure='not_evaluated' (SC-QC-013)"
        )
        if on_scrublet_failure == "raise":
            raise ValueError(message)
        logger.warning(message)
    if "scdblfinder" in detectors:
        table = run_scdblfinder(adata, capture_key, layer, seed, scdblfinder_dbr, rscript, workdir)
        score["scdblfinder"] = table["score"].to_numpy().astype(float)
        call["scdblfinder"] = (table["class"] == "doublet").to_numpy()
        for row in rows.values():
            row["scdblfinder_dbr"] = scdblfinder_dbr
    n_evaluated = np.sum([evaluated[det] for det in detectors], axis=0)
    n_calls = np.sum([call[det] & evaluated[det] for det in detectors], axis=0)
    consensus = (n_calls == n_evaluated) & (n_evaluated >= 2)
    any_call = n_calls >= 1
    for capture, row in rows.items():
        in_capture = captures == capture
        for det in detectors:
            row[f"n_{OBS_PREFIX[det]}"] = int((call[det] & evaluated[det])[in_capture].sum())
        row["n_consensus"] = int(consensus[in_capture].sum())
        row["n_any"] = int(any_call[in_capture].sum())
    for det in detectors:
        adata.obs[f"{OBS_PREFIX[det]}_score"] = score[det]
        adata.obs[f"{OBS_PREFIX[det]}_call"] = call[det] & evaluated[det]
    if "scrublet" in detectors:
        adata.obs["scrublet_evaluated"] = evaluated["scrublet"]
    if "doubletdetection" in detectors:
        adata.obs["dd_label"] = dd_label
    adata.obs["doublet_n_evaluated"] = n_evaluated
    adata.obs["doublet_n_calls"] = n_calls
    adata.obs["doublet_consensus"] = consensus
    adata.obs["doublet_any"] = any_call
    return pl.DataFrame(list(rows.values()))


def combine_calls(*calls: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (consensus, union) of two or more call vectors (1 = doublet); NaN is not called."""
    if len(calls) < 2:
        raise ValueError("combine_calls needs at least two call vectors")
    called = [np.nan_to_num(np.asarray(c, dtype=float), nan=0.0) == 1 for c in calls]
    return np.logical_and.reduce(called), np.logical_or.reduce(called)


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
        Counts without ambient correction, for co-expression (checked with
        ``ambient_method_of_layer`` when ``lineage_markers`` is given).
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
    _warn_if_masking_layer_exists(adata)
    clusters = adata.obs[cluster_key].astype(str).to_numpy()
    captures = adata.obs[capture_key].astype(str).to_numpy()
    any_call = adata.obs["doublet_any"].to_numpy(dtype=bool)
    consensus = adata.obs["doublet_consensus"].to_numpy(dtype=bool)
    capture_rate = {c: float(any_call[captures == c].mean()) for c in np.unique(captures)}
    expected_per_cell = np.array([capture_rate[c] for c in captures])
    n_lineages = None
    if lineage_markers:
        ambient_method_of_layer(adata, layer)
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


def doublet_removal_mask(
    adata: ad.AnnData,
    cluster_key: str,
    review: pl.DataFrame,
    genotype_key: str | None = None,
) -> np.ndarray:
    """Return a mask of cells to remove: ``doublet_consensus``, a ``remove`` cluster, or genotype.

    ``review`` is the output of ``doublet_cluster_review``; set the status of ``review`` clusters
    manually before calling (to ``remove`` or ``keep``) and record the decision.
    ``genotype_key`` names a boolean ``obs`` column with the doublets called by genotype
    demultiplexing in pooled captures (``SC-QC-014``); None for captures that are not pooled.
    """
    if (review["status"] == "review").any():
        logger.warning(
            "clusters with status 'review' remain; they are kept unless set to 'remove'"
        )
    remove_clusters = set(review.filter(pl.col("status") == "remove")["cluster"].to_list())
    in_removed = adata.obs[cluster_key].astype(str).isin(remove_clusters).to_numpy()
    mask = adata.obs["doublet_consensus"].to_numpy(dtype=bool) | in_removed
    if genotype_key is not None:
        mask |= adata.obs[genotype_key].to_numpy(dtype=bool)
    return mask
