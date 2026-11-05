"""Recipe: which layer each step reads, integer counts for count-based DE, protected layers.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §1.1 (``SC-DATA-001``, ``SC-DATA-003``, ``SC-DATA-004``).

``input_layer(adata, step)`` implements the table of ``SC-DATA-004`` (project default; deviations
go to ``PROJECT_CONTEXT.md``):

* ``"qc"`` (cell-level QC metrics) and ``"doublets"`` (doublet detection and co-expression
  review): ``counts``;
* ``"normalization"`` (size factors, normalization, HVG, PCA, clustering, annotation) and
  ``"pseudobulk"`` (pseudobulk DE, abundance, depth matching, empirical power):
  ``counts_ambient_corrected`` if it exists, else ``counts``.

Recipes of the normalization and pseudobulk steps take ``layer=None``, resolved here; an explicit
layer is honoured. QC and doublet recipes keep ``layer="counts"``. The choice is logged.

Count-based DE needs integer counts. SoupX and decontX return non-integer corrected counts;
``require_integer_counts`` rejects them with an error that names the layer and the two remedies:
round explicitly and record it in ``PROJECT_CONTEXT.md`` (e.g. SoupX
``adjustCounts(roundToInt = TRUE)``), or pass ``layer="counts"``. Values are never rounded
silently.

``check_output_layer`` refuses to write a derived layer into a layer of observed counts
(``PROTECTED_LAYERS``, ``SC-DATA-003``) or into the layer being read.
"""

import logging

import anndata as ad
import numpy as np
import scipy.sparse as sp

logger = logging.getLogger(__name__)

COUNTS_LAYER = "counts"
AMBIENT_LAYER = "counts_ambient_corrected"
# Layers of observed counts (SC-DATA-001); derived layers are not written into them (SC-DATA-003).
PROTECTED_LAYERS = (COUNTS_LAYER, AMBIENT_LAYER, "spliced", "unspliced", "ambiguous")
UNCORRECTED_STEPS = ("qc", "doublets")
CORRECTED_STEPS = ("normalization", "pseudobulk")
STEPS = UNCORRECTED_STEPS + CORRECTED_STEPS


def input_layer(adata: ad.AnnData, step: str) -> str:
    """Return the layer that ``step`` reads by default (``SC-DATA-004``) and log the choice.

    Parameters
    ----------
    adata : ad.AnnData
        Data with named layers (``SC-DATA-001``).
    step : str
        One of ``STEPS``: ``"qc"``, ``"doublets"``, ``"normalization"`` or ``"pseudobulk"``.

    Returns
    -------
    str
        ``"counts"`` for QC and doublets; for normalization and pseudobulk,
        ``"counts_ambient_corrected"`` if that layer exists, else ``"counts"``.

    Raises
    ------
    ValueError
        If ``step`` is unknown or the resolved layer is not in ``adata.layers``.
    """
    if step not in STEPS:
        raise ValueError(f"unknown step {step!r}; use one of {STEPS} (SC-DATA-004)")
    if step in CORRECTED_STEPS and AMBIENT_LAYER in adata.layers:
        layer = AMBIENT_LAYER
    else:
        layer = COUNTS_LAYER
    if layer not in adata.layers:
        raise ValueError(
            f"step {step!r} reads layer {layer!r}, which is not in adata.layers "
            "(SC-DATA-001, SC-DATA-004)"
        )
    message = f"step {step!r}: input layer {layer!r} (SC-DATA-004)"
    if step == "pseudobulk" and layer == AMBIENT_LAYER:
        message += f"; check key genes on {COUNTS_LAYER!r}"
    logger.info(message)
    return layer


def require_integer_counts(adata: ad.AnnData, layer: str) -> None:
    """Raise unless ``layers[layer]`` holds non-negative integer counts (input of count-based DE).

    The check is exact (no tolerance): float storage represents integer counts exactly (float32 up
    to 2**24, float64 up to 2**53).

    Raises
    ------
    ValueError
        If the layer is absent, has negative or non-finite values, or non-integer values. The
        message names the layer and the remedies; values are never rounded here.
    """
    if layer not in adata.layers:
        raise ValueError(f"layer {layer!r} not in adata.layers")
    matrix = adata.layers[layer]
    values = np.asarray(matrix.data if sp.issparse(matrix) else matrix)  # native dtype, no copy
    if values.dtype.kind not in "biuf":
        raise ValueError(f"layer {layer!r} has dtype {values.dtype}: not a count matrix")
    if np.any(values < 0) or (values.dtype.kind == "f" and not np.all(np.isfinite(values))):
        raise ValueError(f"layer {layer!r} has negative or non-finite values: not a count matrix")
    if values.dtype.kind == "f" and np.any(values != np.round(values)):
        remedies = (
            "round it explicitly and record it in PROJECT_CONTEXT.md (e.g. SoupX "
            "adjustCounts(roundToInt = TRUE))"
        )
        if layer != COUNTS_LAYER:
            remedies += f", or pass layer={COUNTS_LAYER!r}"
        raise ValueError(
            f"layer {layer!r} has non-integer values; count-based DE needs integers: {remedies} "
            "(SC-DATA-004)"
        )


def check_output_layer(key_added: str, layer: str) -> None:
    """Raise if a derived layer would overwrite observed counts or the layer being read.

    Raises
    ------
    ValueError
        If ``key_added`` is in ``PROTECTED_LAYERS`` or equals ``layer``.
    """
    if key_added in PROTECTED_LAYERS:
        raise ValueError(
            f"key_added={key_added!r} would overwrite a layer of observed counts "
            f"{PROTECTED_LAYERS} (SC-DATA-003); choose another name"
        )
    if key_added == layer:
        raise ValueError(f"key_added={key_added!r} is the input layer; choose another name")
