"""Tests for recipes/scrnaseq/library_qc.py on synthetic STARsolo outputs (nf-core/scrnaseq layout)."""

import logging

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from recipes.scrnaseq.library_qc import (
    FLAG_DIRECTIONS,
    flag_libraries,
    library_qc_table,
    nfcore_scrnaseq_paths,
    read_barcodes_stats,
    read_summary_csv,
    star_version,
    store_library_qc,
)

BARCODE_NAMES = [
    "noNoAdapter",
    "noNoUMI",
    "noNoCB",
    "noNinCB",
    "noNinUMI",
    "noUMIhomopolymer",
    "noNoWLmatch",
    "noTooManyMM",
    "noTooManyWLmatches",
    "yesWLmatchExact",
    "yesOneWLmatchWithMM",
    "yesMultWLmatchWithMM",
]


def summary_text(
    feature="Gene", valid=0.97, in_cells=0.85, cells=5000, saturation=0.6, multi="NoMulti"
):
    f = feature
    rows = [
        ("Number of Reads", 100_000_000),
        ("Reads With Valid Barcodes", valid),
        ("Sequencing Saturation", saturation),
        ("Q30 Bases in CB+UMI", 0.95),
        ("Q30 Bases in RNA read", 0.91),
        ("Reads Mapped to Genome: Unique+Multiple", 0.93),
        ("Reads Mapped to Genome: Unique", 0.88),
        (f"Reads Mapped to {f}: Unique+Multiple {f}", multi),
        (f"Reads Mapped to {f}: Unique {f}", 0.62),
        ("Estimated Number of Cells", cells),
        (f"Unique Reads in Cells Mapped to {f}", 50_000_000),
        ("Fraction of Unique Reads in Cells", in_cells),
        ("Mean Reads per Cell", 10_000),
        ("Median Reads per Cell", 9_000),
        ("UMIs in Cells", 20_000_000),
        ("Mean UMI per Cell", 4_000),
        ("Median UMI per Cell", 3_500),
        (f"Mean {f} per Cell", 1_800),
        (f"Median {f} per Cell", 1_700),
        (f"Total {f} Detected", 25_000),
    ]
    return "".join(f"{name},{value}\n" for name, value in rows)


def barcodes_text(no_match=1_000_000, exact=95_000_000):
    values = dict.fromkeys(BARCODE_NAMES, 1000)
    values["noNoWLmatch"] = no_match
    values["yesWLmatchExact"] = exact
    return "".join(f"{name:>50}{count:>15}\n" for name, count in values.items())


def write_sample(outdir, sample, feature="Gene", **kwargs):
    paths = nfcore_scrnaseq_paths(outdir, sample, feature=feature)
    paths["summary"].parent.mkdir(parents=True)
    barcode_kwargs = {key: kwargs.pop(key) for key in ("no_match", "exact") if key in kwargs}
    paths["summary"].write_text(summary_text(feature=feature, **kwargs), encoding="utf-8")
    paths["barcodes_stats"].write_text(barcodes_text(**barcode_kwargs), encoding="utf-8")
    paths["log_out"].write_text(
        "STAR version=2.7.11b\nSTAR compilation time,server,dir=...\n", encoding="utf-8"
    )
    return paths


def test_nfcore_layout(tmp_path):
    paths = nfcore_scrnaseq_paths(tmp_path, "S1", feature="GeneFull_Ex50pAS")
    assert (
        paths["summary"]
        == tmp_path / "star" / "S1" / "S1.Solo.out" / "GeneFull_Ex50pAS" / "Summary.csv"
    )
    assert paths["barcodes_stats"] == tmp_path / "star" / "S1" / "S1.Solo.out" / "Barcodes.stats"
    assert paths["log_out"] == tmp_path / "star" / "S1" / "S1.Log.out"


def test_summary_detects_feature_and_keeps_nomulti_missing(tmp_path):
    paths = write_sample(tmp_path, "S1", feature="GeneFull_Ex50pAS")
    values = read_summary_csv(paths["summary"])
    assert values["feature"] == "GeneFull_Ex50pAS"
    assert values["median_genes_per_cell"] == 1700
    assert values["frac_mapped_feature_unique"] == pytest.approx(0.62)
    assert values["frac_mapped_feature"] is None  # NoMulti is not 0


def test_summary_rejects_malformed_line(tmp_path):
    path = tmp_path / "Summary.csv"
    path.write_text("Number of Reads 100\n", encoding="utf-8")
    with pytest.raises(ValueError, match="name,value"):
        read_summary_csv(path)


def test_barcodes_stats_and_whitelist_fraction(tmp_path):
    paths = write_sample(tmp_path, "S1", no_match=10_000_000, exact=80_000_000)
    counts = read_barcodes_stats(paths["barcodes_stats"])
    assert counts["noNoWLmatch"] == 10_000_000
    table = library_qc_table({"S1": paths})
    expected = 10_000_000 / (10_000_000 + 80_000_000 + 1000 * 10)
    assert table["frac_cb_no_whitelist_match"][0] == pytest.approx(expected)


def test_star_version_and_missing_log(tmp_path):
    paths = write_sample(tmp_path, "S1")
    assert star_version(paths["log_out"]) == "2.7.11b"
    assert star_version(tmp_path / "absent.Log.out") is None


def test_missing_file_is_an_error(tmp_path):
    paths = write_sample(tmp_path, "S1")
    paths["barcodes_stats"].unlink()
    with pytest.raises(FileNotFoundError, match="barcodes_stats"):
        library_qc_table({"S1": paths})


def five_libraries(tmp_path):
    samples = {
        f"S{i}": write_sample(
            tmp_path,
            f"S{i}",
            valid=0.97 + 0.002 * i,
            in_cells=0.85 + 0.005 * i,
            no_match=1_000_000 + 100_000 * i,
        )
        for i in range(1, 5)
    }
    # S5: wrong whitelist (few valid barcodes, many unmatched) and cell calling picked up empty droplets
    samples["S5"] = write_sample(
        tmp_path,
        "S5",
        valid=0.55,
        in_cells=0.40,
        cells=20_000,
        no_match=40_000_000,
        exact=55_000_000,
    )
    return samples


def test_departing_library_is_flagged(tmp_path, caplog):
    samples = five_libraries(tmp_path)
    loaded = {sample: 10_000 for sample in samples}
    table = library_qc_table(samples, cells_loaded=loaded)
    with caplog.at_level(logging.WARNING):
        flagged, report = flag_libraries(table)
    by_sample = {row["sample_id"]: row for row in flagged.iter_rows(named=True)}
    assert by_sample["S5"]["library_flagged"]
    assert by_sample["S5"]["flag_frac_valid_barcodes"]
    assert by_sample["S5"]["flag_frac_unique_reads_in_cells"]
    assert by_sample["S5"]["flag_frac_cb_no_whitelist_match"]
    assert not any(by_sample[s]["library_flagged"] for s in ["S1", "S2", "S3", "S4"])
    status = dict(zip(report["metric"], report["status"], strict=True))
    assert status["frac_valid_barcodes"] == "evaluated"
    assert "depart from the others" in caplog.text


def test_identical_metric_is_not_evaluated(tmp_path):
    flagged, report = flag_libraries(library_qc_table(five_libraries(tmp_path)))
    status = dict(zip(report["metric"], report["status"], strict=True))
    assert status["q30_rna"] == "not_evaluated_mad_zero"  # same value in every library
    assert flagged["flag_q30_rna"].null_count() == flagged.height


def test_too_few_libraries_and_unknown_loading(tmp_path):
    samples = {s: write_sample(tmp_path, s) for s in ["S1", "S2"]}
    flagged, report = flag_libraries(library_qc_table(samples))
    status = dict(zip(report["metric"], report["status"], strict=True))
    assert status["frac_valid_barcodes"] == "not_evaluated_too_few"
    assert status["cells_over_loaded"] == "not_available"  # no cells_loaded given: missing, not 0
    assert set(status) == set(FLAG_DIRECTIONS)
    assert not flagged["library_flagged"].any()


def test_store_requires_every_sample_and_writes_h5ad(tmp_path):
    samples = five_libraries(tmp_path)
    flagged, _ = flag_libraries(library_qc_table(samples))
    obs = pd.DataFrame(
        {"sample_id": np.repeat(sorted(samples), 3)}, index=[f"c{i}" for i in range(15)]
    )
    adata = ad.AnnData(X=np.ones((15, 2), dtype=np.float32), obs=obs)
    store_library_qc(adata, flagged)
    assert list(adata.uns["library_qc"].index) == sorted(samples)
    # anndata versions still in use refuse nullable string arrays unless this setting is on
    if hasattr(ad.settings, "allow_write_nullable_strings"):
        with ad.settings.override(allow_write_nullable_strings=False):
            adata.write_h5ad(tmp_path / "strict.h5ad")
    adata.write_h5ad(tmp_path / "out.h5ad")
    back = ad.read_h5ad(tmp_path / "out.h5ad")
    assert bool(back.uns["library_qc"].loc["S5", "library_flagged"])

    adata.obs["sample_id"] = ["S9"] + list(adata.obs["sample_id"][1:])
    with pytest.raises(ValueError, match="S9"):
        store_library_qc(adata, flagged)
