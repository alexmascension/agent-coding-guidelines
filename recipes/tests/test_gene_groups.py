"""Tests for recipes/common/gene_groups.py, the shipped lists and their use in qc.py (SC-QC-006)."""

import logging

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pytest
import scipy.sparse as sp

from recipes.common.gene_groups import (
    GROUPS,
    SPECIES,
    excluded_symbols,
    gene_group_path,
    read_gene_group,
)
from recipes.scrnaseq.qc import flag_gene_groups, qc_metrics

ID_PREFIX = {"human": "ENSG", "mouse": "ENSMUSG"}


def test_excluded_symbols_follow_the_guide():
    symbols = [
        "RPL3",
        "RPS6",
        "RPLP0",
        "RPS6KA1",
        "RPS6KB2",
        "MRPL12",
        "MRPS5",
        "HBB",
        "HBA1",
        "HBEGF",
        "HBS1L",
    ]
    assert excluded_symbols(symbols) == ["RPS6KA1", "RPS6KB2", "MRPL12", "MRPS5", "HBEGF", "HBS1L"]
    assert excluded_symbols(["Rps6ka1", "Hbegf", "Rpl3"]) == ["Rps6ka1", "Hbegf"]


def test_read_gene_group_metadata_and_errors(tmp_path):
    path = tmp_path / "g.tsv"
    path.write_text(
        "# species: human\n# source: test\ngene_id\tsymbol\nENSG1\tRPL3\nENSG2\tRPS6\n",
        encoding="utf-8",
    )
    table, meta = read_gene_group(path)
    assert table["gene_id"].to_list() == ["ENSG1", "ENSG2"]
    assert meta == {"species": "human", "source": "test"}
    path.write_text("gene_id\tsymbol\nENSG1\tRPL3\nENSG1\tRPL3\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicated"):
        read_gene_group(path)
    path.write_text("id\tname\nENSG1\tRPL3\n", encoding="utf-8")
    with pytest.raises(ValueError, match="header"):
        read_gene_group(path)
    with pytest.raises(ValueError, match="species"):
        gene_group_path("rat", "ribosomal")


@pytest.mark.parametrize("species", SPECIES)
@pytest.mark.parametrize("group", GROUPS)
def test_shipped_lists(species, group):
    table, meta = read_gene_group(gene_group_path(species, group))
    assert table.height > 0
    assert table["gene_id"].str.starts_with(ID_PREFIX[species]).all()
    assert excluded_symbols(table["symbol"].to_list()) == []
    for key in ("species", "source", "downloaded", "ensembl_release", "hgnc_sha256"):
        assert meta.get(key), f"{key} missing in {species}_{group}"
    assert meta["species"] == species


def _adata():
    var = pd.DataFrame(index=["ENSG_MT1", "ENSG_RPL3", "ENSG_A", "ENSG_B"])
    counts = np.array([[1, 5, 10, 4], [0, 0, 6, 4]], dtype=np.float32)
    adata = ad.AnnData(X=sp.csr_matrix(counts), var=var)
    adata.layers["counts"] = sp.csr_matrix(counts)
    annotation = pl.DataFrame({"gene_id": list(var.index), "chromosome": ["MT", "1", "1", "2"]})
    return adata, annotation


def test_missing_genes_counted_and_empty_group_is_not_zero(caplog):
    adata, annotation = _adata()
    with caplog.at_level(logging.WARNING):
        counts = flag_gene_groups(
            adata,
            annotation,
            ribosomal_ids=["ENSG_RPL3", "ENSG_RPS6"],
            haemoglobin_ids=["ENSG_HBB"],
        )
    assert counts["ribo"] == 1 and counts["ribo_missing"] == 1
    assert counts["hb"] == 0 and counts["hb_missing"] == 1
    assert "no hb gene of the list is in the matrix" in caplog.text
    qc_metrics(adata)
    np.testing.assert_allclose(adata.obs["ribo_fraction"].to_numpy(), [0.25, 0.0])
    assert np.isnan(adata.obs["hb_fraction"]).all()  # not evaluated, never 0


def test_species_uses_shipped_lists():
    adata, annotation = _adata()
    counts = flag_gene_groups(adata, annotation, species="human")
    assert counts["ribo_missing"] > 0  # synthetic IDs: the shipped genes are absent, and counted
