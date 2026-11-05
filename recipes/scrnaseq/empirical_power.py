"""Recipe: empirical power of pseudobulk DE, subsampling donors and cells, held-out reference.

Guidelines: ``SCRNASEQ_GUIDELINES.md`` §2 (``SC-DES-007``) and ``STATISTICS_GUIDELINES.md`` §10
(``STAT-RES-008``), following the design of [Crampton 2026].

Procedure, for one cell type of a large dataset:

1. ``heldout_splits`` divides the units (donors) into a reference and a validation set, stratified
   by the contrast variable, ``n_splits`` times.
2. Reference: pseudobulk DE on all cells of the reference units (units with fewer than
   ``min_cells`` cells left out). The reference set is the genes with ``p_adjusted < alpha``; their
   absolute log2 fold change (shrunken when available) places them in effect-size bins
   (``LFC_BIN_EDGES``); genes below the first edge are in no bin.
3. Designs: for each number of units per group and number of cells per unit (``None`` = all cells
   of each unit), and each of ``n_iterations`` draws, units are drawn per group from the
   validation units with at least that many cells, cells are drawn without replacement within
   each unit, and the same DE is run. A design with too few eligible units is ``infeasible``.
4. Recall per bin: fraction of the reference genes of the bin called (``p_adjusted < alpha``).
   Precision per bin: among the genes called whose own absolute log2 fold change falls in the bin,
   the fraction in the reference set (``p_adjusted < alpha``, any effect size).

Recall and precision against a held-out reference measure agreement with an imperfect reference,
not true power or precision: recall tends to be optimistic (the reference is enriched for large,
well-detected effects) and precision pessimistic (true changes missing from the reference count as
false), but neither is a guaranteed bound. Reference and validation units also share the cohort,
protocol and pipeline, so the estimates describe this cohort. Every design is a full DE run: the
cost is ``n_splits * len(units_per_group) * len(cells_per_unit) * n_iterations`` runs. Seeds derive
from ``seed``, the split, the design and the iteration.

The default layer is resolved with ``recipes/scrnaseq/layers.py`` (``SC-DATA-004``, step
``"pseudobulk"``) and MUST hold integer counts.
"""

import logging

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp

from recipes.scrnaseq.layers import input_layer, require_integer_counts
from recipes.scrnaseq.pseudobulk import ALPHA, MIN_CELLS, pseudobulk_de
from recipes.scrnaseq.units import unit_metadata
from recipes.stats.heldout import REFERENCE_FRACTION, heldout_splits

logger = logging.getLogger(__name__)

N_SPLITS = 3  # held-out splits (Crampton 2026 used 10)
N_ITERATIONS = 5  # draws per design and split (Crampton 2026 used 20)
LFC_BIN_EDGES = (0.2, 1.0, 2.0, 3.0)  # absolute log2 fold change; the last bin is open
RUN_SCHEMA = {
    "split": pl.Int64,
    "units_per_group": pl.Int64,
    "cells_per_unit": pl.Int64,
    "iteration": pl.Int64,
    "lfc_bin": pl.Utf8,
    "status": pl.Utf8,
    "n_reference": pl.Int64,
    "n_called": pl.Int64,
    "n_recovered": pl.Int64,
    "recall": pl.Float64,
    "precision": pl.Float64,
}


def _bin_labels(edges: tuple[float, ...]) -> list[str]:
    """Labels of the effect-size bins, e.g. ``"0.2-1"`` and ``">3"``."""
    fmt = [f"{e:g}" for e in edges]
    return [f"{a}-{b}" for a, b in zip(fmt[:-1], fmt[1:], strict=True)] + [f">{fmt[-1]}"]


def _assign_bins(abs_lfc: np.ndarray, edges: tuple[float, ...]) -> np.ndarray:
    """Bin label of each absolute log fold change; None below the first edge or not finite."""
    labels = np.array(_bin_labels(edges), dtype=object)
    which = np.searchsorted(np.asarray(edges), abs_lfc, side="right") - 1
    out = np.full(abs_lfc.size, None, dtype=object)
    ok = np.isfinite(abs_lfc) & (which >= 0)
    out[ok] = labels[which[ok]]
    return out


def _pseudobulk(
    counts: sp.csr_matrix,
    cells_by_unit: dict[str, np.ndarray],
    meta: pd.DataFrame,
    genes: pd.Index,
    cell_type: str,
) -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    """Sum the given cells of each unit into the input of ``pseudobulk_de``."""
    names = list(cells_by_unit)
    summed = np.vstack([np.asarray(counts[cells_by_unit[u]].sum(axis=0)).ravel() for u in names])
    frame = pd.DataFrame(summed, index=names, columns=genes)
    unit_meta = meta.loc[names].copy()
    unit_meta["n_cells"] = [cells_by_unit[u].size for u in names]
    return {cell_type: (frame, unit_meta)}


def empirical_power(
    adata: ad.AnnData,
    unit_key: str,
    celltype_key: str,
    cell_type: str,
    covariates: list[str],
    design: str,
    contrast: tuple[str, str, str],
    units_per_group: list[int],
    cells_per_unit: list[int | None],
    n_splits: int = N_SPLITS,
    n_iterations: int = N_ITERATIONS,
    reference_fraction: float = REFERENCE_FRACTION,
    lfc_bin_edges: tuple[float, ...] = LFC_BIN_EDGES,
    alpha: float = ALPHA,
    min_cells: int = MIN_CELLS,
    layer: str | None = None,
    seed: int = 0,
    n_cpus: int = 1,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Recall and precision of pseudobulk DE over a grid of designs, against a held-out reference.

    Parameters
    ----------
    adata : ad.AnnData
        Cells after QC; integer counts in ``layer``.
    unit_key : str
        ``obs`` column with the experimental unit (donor).
    celltype_key : str
        ``obs`` column with the annotation.
    cell_type : str
        Cell type analysed.
    covariates : list[str]
        Unit-level ``obs`` columns of the design (the contrast variable included).
    design : str
        Design formula of the real analysis.
    contrast : tuple[str, str, str]
        ``(variable, tested_level, reference_level)``; the splits are stratified by ``variable``.
    units_per_group : list[int]
        Numbers of units drawn per compared group.
    cells_per_unit : list[int | None]
        Numbers of cells drawn per unit; ``None`` uses all cells of each unit.
    n_splits : int
        Held-out splits.
    n_iterations : int
        Draws per design and split.
    reference_fraction : float
        Share of units in the reference.
    lfc_bin_edges : tuple[float, ...]
        Edges of the effect-size bins (absolute log2 fold change).
    alpha : float
        FDR level of every DE call.
    min_cells : int
        Minimum cells for a unit to enter the reference DE.
    layer : str | None
        Integer counts layer; None resolves it with ``input_layer(adata, "pseudobulk")``.
    seed : int
        Base seed.
    n_cpus : int
        CPUs for pydeseq2.

    Returns
    -------
    runs : pl.DataFrame
        One row per split, design, iteration and bin (``RUN_SCHEMA``); ``cells_per_unit`` is
        missing for all cells. ``status`` is ``ok``, ``infeasible`` (too few eligible units),
        ``not_testable`` (the DE of the design could not run) or ``reference_not_testable`` (the
        reference DE of the split could not run: no reference unit with ``min_cells`` cells, too
        few units per group, or no gene passing the filter; the split is not used).
    summary : pl.DataFrame
        Per design and bin: ``mean_recall``, ``mean_precision``, ``n_runs`` (``ok`` runs),
        ``n_infeasible`` and ``n_reference_not_testable``.

    Raises
    ------
    ValueError
        If ``cell_type`` is absent, a grid value is below 2, or the contrast levels are missing.
    """
    variable, tested, reference_level = contrast
    if any(n < 2 for n in units_per_group):
        raise ValueError(f"units_per_group values must be at least 2, got {units_per_group}")
    if any(c is not None and c < 1 for c in cells_per_unit):
        raise ValueError(f"cells_per_unit values must be at least 1 or None, got {cells_per_unit}")
    in_type = (adata.obs[celltype_key].astype(str) == cell_type).to_numpy()
    if not in_type.any():
        raise ValueError(f"cell type {cell_type!r} not in {celltype_key!r}")
    obs = adata.obs.loc[in_type]
    meta = unit_metadata(obs, unit_key, covariates)
    levels = meta[variable].astype(str)
    if not {tested, reference_level} <= set(levels):
        raise ValueError(f"levels {tested!r} and {reference_level!r} must be in {variable!r}")
    layer = input_layer(adata, "pseudobulk") if layer is None else layer
    require_integer_counts(adata, layer)
    layer_counts = adata.layers[layer]
    layer_counts = (
        layer_counts.tocsr() if sp.issparse(layer_counts) else sp.csr_matrix(layer_counts)
    )
    counts = layer_counts[np.flatnonzero(in_type)]
    cell_units = obs[unit_key].astype(str).to_numpy()
    cells_of = {u: np.flatnonzero(cell_units == u) for u in meta.index}
    bins = _bin_labels(lfc_bin_edges)
    rows = []
    splits = heldout_splits(
        meta.index.to_numpy(), levels.to_numpy(), n_splits, reference_fraction, seed
    )
    for s, (reference_units, validation_units) in enumerate(splits):
        ref_units = [u for u in reference_units if cells_of[u].size >= min_cells]
        ref_results, ref_not_testable = pl.DataFrame(), pl.DataFrame()
        if ref_units:
            ref_pb = _pseudobulk(
                counts, {u: cells_of[u] for u in ref_units}, meta, adata.var_names, cell_type
            )
            ref_results, ref_not_testable = pseudobulk_de(
                ref_pb, design, contrast, alpha=alpha, n_cpus=n_cpus
            )
        if not ref_units or ref_not_testable.height or ref_results.height == 0:
            reason = (
                ref_not_testable["reason"][0]
                if ref_not_testable.height
                else f"{len(ref_units)} reference units with at least {min_cells} cells"
            )
            logger.warning(f"split {s}: reference DE not testable ({reason}); split not used")
            rows += [
                {
                    "split": s,
                    "units_per_group": n_units,
                    "cells_per_unit": n_cells,
                    "iteration": it,
                    "lfc_bin": b,
                    "status": "reference_not_testable",
                }
                for n_units in units_per_group
                for n_cells in cells_per_unit
                for it in range(n_iterations)
                for b in bins
            ]
            continue
        ref = ref_results.to_pandas()
        ref = ref[ref["p_adjusted"] < alpha]
        lfc = ref["log2_fold_change_shrunk"].where(
            np.isfinite(ref["log2_fold_change_shrunk"]), ref["log2_fold_change"]
        )
        ref_bins = _assign_bins(np.abs(lfc.to_numpy(dtype=float)), lfc_bin_edges)
        reference_set = set(ref["gene_id"])
        reference_by_bin = {b: set(ref["gene_id"][ref_bins == b]) for b in bins}
        logger.info(
            f"split {s}: reference of {len(ref_units)} units, {len(reference_set)} genes called"
        )
        if not reference_set:
            logger.warning(
                f"split {s}: the reference calls no gene; recall is undefined and every call of "
                "a design counts against precision"
            )
        val_levels = levels.loc[validation_units]
        for d_i, n_units in enumerate(units_per_group):
            for c_i, n_cells in enumerate(cells_per_unit):
                need = 1 if n_cells is None else n_cells
                eligible = {
                    lv: [u for u in val_levels.index[val_levels == lv] if cells_of[u].size >= need]
                    for lv in (tested, reference_level)
                }
                feasible = all(len(v) >= n_units for v in eligible.values())
                for it in range(n_iterations):
                    base = {
                        "split": s,
                        "units_per_group": n_units,
                        "cells_per_unit": n_cells,
                        "iteration": it,
                    }
                    if not feasible:
                        rows += [{**base, "lfc_bin": b, "status": "infeasible"} for b in bins]
                        continue
                    rng = np.random.default_rng([seed, s, d_i, c_i, it])
                    drawn = {}
                    for lv in (tested, reference_level):
                        for u in rng.choice(eligible[lv], size=n_units, replace=False):
                            cells = cells_of[u]
                            drawn[u] = (
                                cells
                                if n_cells is None
                                else np.sort(rng.choice(cells, size=n_cells, replace=False))
                            )
                    pb = _pseudobulk(counts, drawn, meta, adata.var_names, cell_type)
                    results, not_testable = pseudobulk_de(
                        pb, design, contrast, alpha=alpha, n_cpus=n_cpus, shrink=False
                    )
                    if not_testable.height or results.height == 0:
                        rows += [{**base, "lfc_bin": b, "status": "not_testable"} for b in bins]
                        continue
                    res = results.to_pandas()
                    called = res[res["p_adjusted"] < alpha]
                    called_bins = _assign_bins(
                        np.abs(called["log2_fold_change"].to_numpy(dtype=float)), lfc_bin_edges
                    )
                    called_set = set(called["gene_id"])
                    for b in bins:
                        truth = reference_by_bin[b]
                        in_bin = set(called["gene_id"][called_bins == b])
                        recovered = len(truth & called_set)
                        rows.append(
                            {
                                **base,
                                "lfc_bin": b,
                                "status": "ok",
                                "n_reference": len(truth),
                                "n_called": len(in_bin),
                                "n_recovered": recovered,
                                "recall": recovered / len(truth) if truth else None,
                                "precision": (
                                    len(in_bin & reference_set) / len(in_bin) if in_bin else None
                                ),
                            }
                        )
    runs = pl.DataFrame([{c: r.get(c) for c in RUN_SCHEMA} for r in rows], schema=RUN_SCHEMA)
    keys = ["units_per_group", "cells_per_unit", "lfc_bin"]
    summary = (
        runs.group_by(keys, maintain_order=True)
        .agg(
            mean_recall=pl.col("recall").filter(pl.col("status") == "ok").mean(),
            mean_precision=pl.col("precision").filter(pl.col("status") == "ok").mean(),
            n_runs=(pl.col("status") == "ok").sum(),
            n_infeasible=(pl.col("status") == "infeasible").sum(),
            n_reference_not_testable=(pl.col("status") == "reference_not_testable").sum(),
        )
        .sort(keys, nulls_last=True)
    )
    return runs, summary
