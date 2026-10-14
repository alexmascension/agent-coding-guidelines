"""Recipe: read gene-level tables (nf-core/rnaseq or ``tximport_lengthscaled.R``) into samples x genes matrices.

Guideline: ``RNASEQ_BULK_GUIDELINES.md`` §4 (quantification): for pydeseq2 use the ``lengthScaledTPM`` counts
(nf-core: ``salmon.merged.gene_counts_length_scaled.tsv``), never TPM.

``read_gene_table`` returns samples x genes with stable Ensembl IDs (versions stripped) and, separately, the
``gene_name`` column when present; ``to_integer_counts`` rounds length-scaled counts (non-integer by construction) as
DESeq2 does with ``DESeqDataSetFromMatrix(round(...))``.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from recipes.common.gene_ids import strip_version

ANNOTATION_COLUMNS = ("gene_id", "gene_name")


def read_gene_table(path: str | Path) -> tuple[pd.DataFrame, pd.Series | None]:
    """Read a genes x samples TSV and return (samples x genes matrix, gene_name per gene or None).

    Raises
    ------
    ValueError
        If ``gene_id`` is missing or IDs are not unique after removing versions.
    """
    table = pd.read_csv(path, sep="\t")
    if "gene_id" not in table.columns:
        raise ValueError(f"{path}: no gene_id column")
    table["gene_id"] = table["gene_id"].map(strip_version)
    if table["gene_id"].duplicated().any():
        raise ValueError(f"{path}: gene IDs not unique after removing versions")
    table = table.set_index("gene_id")
    names = table["gene_name"] if "gene_name" in table.columns else None
    values = table.drop(columns=[c for c in ANNOTATION_COLUMNS if c in table.columns])
    return values.T.astype(float), names


def to_integer_counts(counts: pd.DataFrame) -> pd.DataFrame:
    """Round length-scaled counts to integers; refuses negative or non-finite values."""
    values = counts.to_numpy()
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("counts contain negative or non-finite values")
    return counts.round().astype(np.int64)


def mean_gene_length(length_table: pd.DataFrame) -> pd.Series:
    """Average effective gene length across samples (samples x genes input), for length-bias checks."""
    return length_table.mean(axis=0).rename("length")
