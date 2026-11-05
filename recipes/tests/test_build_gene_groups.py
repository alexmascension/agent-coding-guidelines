"""Tests of recipes/common/build_gene_groups.py with responses written by hand (no network)."""

import hashlib
import json

import pytest

from recipes.common import build_gene_groups as bgg

HGNC = {
    "response": {
        "docs": [
            {
                "symbol": "RPL3",
                "status": "Approved",
                "locus_group": "protein-coding gene",
                "ensembl_gene_id": "ENSG00000100316",
                "gene_group": ["L ribosomal proteins"],
                "mgd_id": ["MGI:1"],
            },
            {
                "symbol": "RPS6",
                "status": "Approved",
                "locus_group": "protein-coding gene",
                "ensembl_gene_id": "ENSG00000137154",
                "gene_group": ["S ribosomal proteins"],
                "mgd_id": ["MGI:2"],
            },
            {
                "symbol": "HBB",
                "status": "Approved",
                "locus_group": "protein-coding gene",
                "ensembl_gene_id": "ENSG00000244734",
                "gene_group": ["Hemoglobin subunits"],
                "mgd_id": ["MGI:3"],
            },
            {
                "symbol": "RPL3P1",
                "status": "Entry Withdrawn",
                "locus_group": "protein-coding gene",
                "ensembl_gene_id": "ENSG00000000001",
                "gene_group": ["L ribosomal proteins"],
                "mgd_id": [],
            },
        ]
    }
}
BIOMART = {
    "MGI:1,MGI:2": "Gene stable ID\tGene name\tGene type\nENSMUSG01\tRpl3\tprotein_coding\n"
    "ENSMUSG02\tRps6\tprotein_coding\nENSMUSG09\tRps6-ps1\tprocessed_pseudogene\n",
    "MGI:3": "Gene stable ID\tGene name\tGene type\nENSMUSG03\tHbb-bs\tprotein_coding\n",
}


class _Response:
    def __init__(self, content: bytes):
        self.content = content


def _fake_get(calls: list):
    def get(url, **kwargs):
        calls.append(url)
        if url == bgg.HGNC_URL:
            return _Response(json.dumps(HGNC).encode())
        if url == f"{bgg.HGNC_REST}/info":
            return _Response(b'{"lastModified": "2026-11-01T00:00:00Z"}')
        if url == f"{bgg.ENSEMBL_REST}/info/data":
            return _Response(b'{"releases": [115]}')
        if url == bgg.BIOMART:
            mgi = kwargs["params"]["query"].split('name="mgi_id" value="')[1].split('"')[0]
            return _Response(BIOMART[mgi].encode())
        raise AssertionError(url)

    return get


def _no_network(url, **kwargs):
    raise AssertionError(f"network request to {url}")


@pytest.fixture
def built(tmp_path, monkeypatch):
    """One build into ``out`` with its snapshot in ``snap``."""
    calls = []
    monkeypatch.setattr(bgg, "_get", _fake_get(calls))
    out, snap = tmp_path / "out", tmp_path / "snap"
    assert bgg.main(["--out", str(out), "--snapshot-dir", str(snap)]) == 0
    monkeypatch.setattr(bgg, "_get", _no_network)
    return {"out": out, "snap": snap, "calls": calls, "tmp": tmp_path}


def _metadata(path):
    return dict(
        line[2:].split(": ", 1) for line in path.read_text().splitlines() if line.startswith("# ")
    )


def test_build_writes_lists_and_a_checksummed_snapshot(built):
    out, snap = built["out"], built["snap"]
    assert sorted(p.name for p in out.iterdir()) == sorted(bgg.OUTPUT_NAMES)
    assert (out / "human_ribosomal.tsv").read_text().splitlines()[-2:] == [
        "ENSG00000100316\tRPL3",
        "ENSG00000137154\tRPS6",
    ]
    assert (out / "mouse_ribosomal.tsv").read_text().splitlines()[-1] == "ENSMUSG02\tRps6"
    index = json.loads((snap / bgg.SNAPSHOT_INDEX).read_text())
    assert sorted(index["files"]) == [
        "biomart_haemoglobin.tsv",
        "biomart_ribosomal.tsv",
        "ensembl_info_data.json",
        "hgnc_info.json",
        "hgnc_protein_coding.json",
    ]
    for entry in index["files"].values():
        assert hashlib.sha256((snap / entry["file"]).read_bytes()).hexdigest() == entry["sha256"]
    meta = _metadata(out / "human_ribosomal.tsv")
    raw = (snap / "hgnc_protein_coding.json").read_bytes()
    assert meta["hgnc_sha256"] == hashlib.sha256(raw).hexdigest()
    sums = (snap / bgg.SNAPSHOT_SUMS).read_bytes()
    assert meta["sources_sha256"] == hashlib.sha256(sums).hexdigest()
    assert meta["ensembl_release"] == "115" and meta["downloaded"] == index["downloaded"]
    assert not list(out.glob(".*.tmp")) and not list(snap.glob(".*.tmp"))


def test_existing_lists_are_replaced_only_with_overwrite(built, monkeypatch):
    out, snap = built["out"], built["snap"]
    before = {p.name: p.read_bytes() for p in out.iterdir()}
    assert bgg.main(["--out", str(out), "--snapshot-dir", str(built["tmp"] / "s2")]) == 1
    assert {p.name: p.read_bytes() for p in out.iterdir()} == before  # nothing fetched or written
    monkeypatch.setattr(bgg, "_get", _fake_get([]))
    assert bgg.main(["--out", str(built["tmp"] / "o2"), "--snapshot-dir", str(snap)]) == 1
    assert bgg.main(["--out", str(out), "--snapshot-dir", str(snap), "--overwrite"]) == 0


def test_rebuild_from_snapshot_is_offline_and_identical(built):
    out2 = built["tmp"] / "rebuilt"
    assert bgg.main(["--from-snapshot", str(built["snap"]), "--out", str(out2)]) == 0
    for name in bgg.OUTPUT_NAMES:
        assert (out2 / name).read_bytes() == (built["out"] / name).read_bytes()


def test_a_changed_snapshot_file_is_refused(built):
    (built["snap"] / "biomart_ribosomal.tsv").write_bytes(
        b"Gene stable ID\tGene name\tGene type\n"
    )
    with pytest.raises(ValueError, match="SHA-256 differs"):
        bgg.main(["--from-snapshot", str(built["snap"]), "--out", str(built["tmp"] / "o3")])
    assert not (built["tmp"] / "o3").exists()


def test_failed_write_keeps_the_previous_file(tmp_path, monkeypatch):
    path = tmp_path / "list.tsv"
    path.write_bytes(b"old\n")

    def fail(src, dst):
        raise OSError("disk full")

    monkeypatch.setattr(bgg.os, "replace", fail)
    with pytest.raises(OSError, match="disk full"):
        bgg.write_atomic(path, b"new\n")
    assert path.read_bytes() == b"old\n"
    assert [p.name for p in tmp_path.iterdir()] == ["list.tsv"]
