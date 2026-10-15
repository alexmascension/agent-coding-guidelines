r"""Recipe: pooled size factors (scranPY; R scran check), median factors, shifted-log transform.

Guideline: ``SCRNASEQ_GUIDELINES.md`` §5 (normalization and feature selection).

Procedure
---------
* ``scranpy_size_factors`` (default): pooled deconvolution size factors [Lun 2016] with ``scranPY``
  (``scranPY.compute_sum_factors``, a Python port of ``scran::computeSumFactors``) [Fortmann 2023].
  Clusters for
  pooling come from ``cluster_key`` (e.g. a first Leiden clustering); without clusters all cells
  are pooled together.
* ``scran_size_factors``: the same through R (``scran_size_factors.R``: ``quickCluster`` +
  ``calculateSumFactors``).
  Use it to cross-check scranPY on a subset when a project adopts it, or when the dense matrix does
  not fit in memory.
* ``median_size_factors``: $s_i = \text{total}_i / \text{median}(\text{total})$; with
  ``log_normalize`` this is
  numerically equal to ``scanpy.pp.normalize_total(target_sum=None)`` followed by ``log1p`` (each
  cell scaled to the
  median total).
* ``log_normalize``: ``layers["log1p_norm"] = log1p(counts / s_i)`` with the factors **as given**
  (no re-centring).
  scran and scranPY factors are returned centred to mean 1 (the scran convention); median factors
  are centred at the
  median cell. The centring sets the effective pseudocount, so it is part of the method and is
  recorded.
* All factors are validated: one per cell, finite and positive.

Installation and limits of scranPY
----------------------------------
* Install from GitHub pinned to a commit: ``pip install
  git+https://github.com/sfortma2/scranPY.git@<commit>``
  (``SCRANPY_COMMIT``). Its distribution name normalizes to ``scranpy`` and collides with the
  unrelated Bioconductor
  ``scranpy`` package (libscran bindings): both cannot be installed in the same environment.
* It needs a dense ``adata.X``; memory is ``n_cells * n_genes * 8`` bytes. ``scranpy_size_factors``
  refuses above
  ``MAX_DENSE_GB`` instead of failing mid-run.
* On simulated data with composition differences between cell types (2,400 cells, 2,000 genes) its
  factors correlated
  at 0.995 (log scale) with R scran, and both were closer to the true factors than library size.
"""

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import polars as pl
import scipy.io
import scipy.sparse as sp

from recipes.common.resources import package_script

logger = logging.getLogger(__name__)

SCRAN_MIN_MEAN = 0.1  # scran recommendation for UMI counts
SCRANPY_COMMIT = "62292a039158fce0123bdec9a5e6b73ecf5db123"  # sfortma2/scranPY, version 0.1
MAX_DENSE_GB = 32.0  # project-adjustable memory guard for the dense copy scranPY needs


def _totals(adata: ad.AnnData, layer: str) -> np.ndarray:
    return np.asarray(adata.layers[layer].sum(axis=1)).ravel().astype(float)


def validate_size_factors(
    size_factors: np.ndarray, n_cells: int, source: str = "input"
) -> np.ndarray:
    """Return factors as float array after checking one finite, positive value per cell.

    Raises
    ------
    ValueError
        If the length differs from ``n_cells`` or any factor is NaN, infinite, zero or negative.
    """
    factors = np.asarray(size_factors, dtype=float).ravel()
    if factors.shape != (n_cells,):
        raise ValueError(f"{source}: {factors.size} size factors for {n_cells} cells")
    bad = ~np.isfinite(factors) | (factors <= 0)
    if bad.any():
        raise ValueError(f"{source}: {int(bad.sum())} size factors are not finite and positive")
    return factors


def median_size_factors(adata: ad.AnnData, layer: str = "counts") -> np.ndarray:
    """Size factors proportional to library size, scaled so that the median cell has factor 1.

    Raises
    ------
    ValueError
        If any cell has zero counts.
    """
    totals = _totals(adata, layer)
    if not np.all(np.isfinite(totals)) or np.any(totals <= 0):
        raise ValueError(
            f"{int(np.sum(~(totals > 0)))} cells with zero or invalid totals; filter them first"
        )
    return totals / np.median(totals)


def scranpy_size_factors(
    adata: ad.AnnData,
    layer: str = "counts",
    cluster_key: str | None = None,
    min_mean: float = SCRAN_MIN_MEAN,
    max_dense_gb: float = MAX_DENSE_GB,
    parallelize: bool = False,
) -> np.ndarray:
    """Pooled size factors with scranPY, centred to mean 1.

    Parameters
    ----------
    adata : ad.AnnData
        Cells after QC; raw counts in ``layer``. Not modified.
    layer : str
        Layer with raw counts.
    cluster_key : str | None
        ``obs`` column with clusters for pooling (recommended when cell types differ in RNA
        content).
    min_mean : float
        Minimum mean count of genes used for pooling (``min.mean`` in scran).
    max_dense_gb : float
        Refuse to build a dense matrix larger than this.
    parallelize : bool
        scranPY multiprocessing over clusters.

    Returns
    -------
    np.ndarray
        Size factor per cell, in the order of ``adata.obs_names``.

    Raises
    ------
    MemoryError
        If the dense copy would exceed ``max_dense_gb``.
    ValueError
        If scranPY returns non-positive factors.
    """
    import contextlib
    import io

    import scranPY

    dense_gb = adata.n_obs * adata.n_vars * 8 / 1e9
    if dense_gb > max_dense_gb:
        raise MemoryError(
            f"dense copy needs {dense_gb:.1f} GB (> {max_dense_gb}); use scran_size_factors (R)"
        )
    counts = adata.layers[layer]
    dense = counts.toarray() if sp.issparse(counts) else np.asarray(counts)
    obs = adata.obs[[cluster_key]].copy() if cluster_key else None
    tmp = ad.AnnData(X=dense.astype(np.float64), obs=obs)
    tmp.obs_names = adata.obs_names
    with contextlib.redirect_stdout(
        io.StringIO()
    ) as captured:  # scranPY prints progress with print()
        factors = scranPY.compute_sum_factors(
            tmp,
            clusters=cluster_key,
            min_mean=min_mean,
            parallelize=parallelize,
            algorithm="CVXPY",
            plotting=False,
            stopwatch=False,
            normalize_counts=False,
        )
    logger.debug(f"scranPY output: {captured.getvalue()}")
    factors = validate_size_factors(factors, adata.n_obs, source="scranPY")
    return factors / factors.mean()


def scran_size_factors(
    adata: ad.AnnData,
    layer: str = "counts",
    cluster_key: str | None = None,
    min_mean: float = SCRAN_MIN_MEAN,
    seed: int = 0,
    rscript: str = "Rscript",
    workdir: str | Path | None = None,
) -> pl.DataFrame:
    """Pooled size factors from scran through ``scran_size_factors.R``.

    Parameters
    ----------
    adata : ad.AnnData
        Cells after QC; raw integer counts in ``layer``.
    layer : str
        Layer with raw counts.
    cluster_key : str | None
        ``obs`` column with clusters for pooling; None lets ``quickCluster`` define them.
    min_mean : float
        ``min.mean`` of ``calculateSumFactors``.
    seed : int
        Seed for ``quickCluster``.
    rscript : str
        ``Rscript`` executable of an environment with scran, Matrix and optparse.
    workdir : str | Path | None
        Parent directory for the exchange files (a unique subdirectory is created and kept, with
        the R provenance
        file); None uses a temporary directory deleted afterwards.

    Returns
    -------
    pl.DataFrame
        ``barcode``, ``size_factor`` (mean 1), ``cluster``, in the order of ``adata.obs_names``.

    Raises
    ------
    RuntimeError
        If the R script fails (its stderr is included).
    """
    if shutil.which(rscript) is None:
        raise RuntimeError(f"{rscript} not found; run in an environment with R and scran")
    if workdir is not None:
        Path(workdir).mkdir(parents=True, exist_ok=True)
    tmp = Path(
        tempfile.mkdtemp(prefix="scran_", dir=workdir)
    )  # unique per run; never reuses earlier files
    r_script = package_script("recipes.scrnaseq", "scran_size_factors.R")
    try:
        counts = adata.layers[layer]
        counts = counts if sp.issparse(counts) else sp.csr_matrix(counts)
        scipy.io.mmwrite(tmp / "counts.mtx", counts.T.tocoo())
        (tmp / "barcodes.txt").write_text("\n".join(adata.obs_names) + "\n")
        command = [
            rscript,
            str(r_script),
            "--counts",
            str(tmp / "counts.mtx"),
            "--barcodes",
            str(tmp / "barcodes.txt"),
            "--min-mean",
            str(min_mean),
            "--seed",
            str(seed),
            "--out",
            str(tmp / "size_factors.tsv"),
        ]
        if cluster_key is not None:
            (tmp / "clusters.txt").write_text("\n".join(adata.obs[cluster_key].astype(str)) + "\n")
            command += ["--clusters", str(tmp / "clusters.txt")]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"scran_size_factors.R failed:\n{result.stderr}")
        table = pl.read_csv(
            tmp / "size_factors.tsv", separator="\t", schema_overrides={"cluster": pl.Utf8}
        )
    finally:
        if workdir is None:
            shutil.rmtree(tmp, ignore_errors=True)
    if table["barcode"].to_list() != list(adata.obs_names):
        raise RuntimeError("barcode order returned by R differs from adata.obs_names")
    validate_size_factors(table["size_factor"].to_numpy(), adata.n_obs, source="scran")
    return table


def log_normalize(
    adata: ad.AnnData,
    size_factors: np.ndarray,
    layer: str = "counts",
    key_added: str = "log1p_norm",
) -> None:
    """Store ``log1p(counts / s)`` in ``layers[key_added]`` and ``obs["size_factor"]`` (in place).

    Raises
    ------
    ValueError
        If factors are not positive or their length differs from the number of cells.
    """
    factors = validate_size_factors(size_factors, adata.n_obs)
    counts = adata.layers[layer]
    if sp.issparse(counts):
        scaled = sp.diags(1.0 / factors) @ counts.tocsr().astype(np.float64)
        scaled.data = np.log1p(scaled.data)
    else:
        scaled = np.log1p(np.asarray(counts, dtype=float) / factors[:, None])
    adata.layers[key_added] = scaled
    adata.obs["size_factor"] = factors
