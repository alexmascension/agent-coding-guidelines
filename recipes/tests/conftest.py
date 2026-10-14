"""Synthetic data shared by the recipe tests (no real data, no network)."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

GTF_LINES = [
    "#!genome-build GRCh38",
    '1\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000000001"; gene_version "3"; gene_name "GFAP"; '
    'gene_biotype "protein_coding";',
    '1\tensembl\tgene\t200\t300\t.\t+\t.\tgene_id "ENSG00000000002"; gene_version "1"; gene_name "AQP4"; '
    'gene_biotype "protein_coding";',
    '1\tensembl\ttranscript\t200\t300\t.\t+\t.\tgene_id "ENSG00000000002"; transcript_id "ENST1";',
    '2\tensembl\tgene\t1\t100\t.\t-\t.\tgene_id "ENSG00000000003"; gene_version "7"; gene_name "MTOR"; '
    'gene_biotype "protein_coding";',
    'MT\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000000004"; gene_version "1"; gene_name "MT-CO1"; '
    'gene_biotype "protein_coding";',
    'X\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000000005"; gene_version "2"; gene_name "DUP"; '
    'gene_biotype "protein_coding";',
    'Y\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000000006"; gene_version "2"; gene_name "DUP"; '
    'gene_biotype "protein_coding";',
    '3\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000000007"; gene_version "1"; gene_biotype "lncRNA";',
]


@pytest.fixture
def gtf_path(tmp_path):
    """Small Ensembl-style GTF with a mitochondrial gene, a duplicated symbol and a gene without symbol."""
    path = tmp_path / "toy.gtf"
    path.write_text("\n".join(GTF_LINES) + "\n")
    return path


def simulate_counts(seed: int = 0, n_genes: int = 300, cells_per_type: int = 60, de_gene_fold: float = 4.0):
    """Negative-binomial counts: 3 cell types x 6 samples (3 control, 3 case), with markers and one DE gene.

    Genes 0-9 mark type A, 10-19 type B, 20-29 type C; genes 30-34 are mitochondrial; gene 40 is up-regulated
    ``de_gene_fold`` times in type A of case samples. Spliced/unspliced layers split
    each count binomially (30% unspliced).
    """
    rng = np.random.default_rng(seed)
    types, samples = ["A", "B", "C"], [f"s{i}" for i in range(6)]
    condition = {s: ("case" if i >= 3 else "control") for i, s in enumerate(samples)}
    base = rng.gamma(2.0, 1.0, size=n_genes)
    rows, obs = [], []
    for t_index, cell_type in enumerate(types):
        mean = base.copy()
        mean[t_index * 10 : (t_index + 1) * 10] *= 15
        mean[30:35] = 3.0
        for sample in samples:
            sample_mean = mean * rng.lognormal(0, 0.1, size=n_genes)
            if cell_type == "A" and condition[sample] == "case":
                sample_mean[40] *= de_gene_fold
            n = cells_per_type // len(samples) + 5
            depth = rng.lognormal(0, 0.3, size=(n, 1))
            mu = sample_mean[None, :] * depth
            rows.append(rng.negative_binomial(5, 5 / (5 + mu)))
            obs += [
                {"cell_type": cell_type, "sample": sample, "condition": condition[sample], "capture": sample}
                for _ in range(n)
            ]
    counts = np.vstack(rows).astype(np.float32)
    obs = pd.DataFrame(obs, index=[f"cell{i}" for i in range(counts.shape[0])])
    var = pd.DataFrame(index=[f"ENSG{i:011d}.1" for i in range(n_genes)])
    adata = ad.AnnData(X=sp.csr_matrix(counts), obs=obs, var=var)
    adata.layers["counts"] = sp.csr_matrix(counts)
    unspliced = rng.binomial(counts.astype(np.int64), 0.3).astype(np.float32)
    adata.layers["unspliced"] = sp.csr_matrix(unspliced)
    adata.layers["spliced"] = sp.csr_matrix(counts - unspliced)
    return adata


@pytest.fixture
def sim_adata():
    """Simulated count data (see ``simulate_counts``)."""
    return simulate_counts()
