"""Smoke and unit tests for the panel builder.

These don't run the full 3-hour stream — they validate the local glue (VCF
parsing, dosage decoding, marker-allele filters, write_compact_panel) using
small synthetic inputs.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomi_ancestry_panel import builder, naming


def test_dosages_from_samples_handles_phased_and_unphased() -> None:
    fmt = "GT:DP"
    fields = ["0|0:30", "0/1:25", "1|1:40", "0|.:10", "./.:8", "1/0:20"]
    out = builder._dosages_from_samples(fmt, fields)
    assert out == [0.0, 1.0, 2.0, None, None, 1.0]


def test_dosages_from_samples_returns_all_none_when_no_gt() -> None:
    out = builder._dosages_from_samples("DP:GQ", ["30:99", "25:80"])
    assert out == [None, None]


def test_usable_marker_alleles_rejects_indels_and_palindromes() -> None:
    assert builder._usable_marker_alleles("A", "G")
    assert builder._usable_marker_alleles("C", "T")
    assert not builder._usable_marker_alleles("AC", "A")        # indel
    assert not builder._usable_marker_alleles("A", "AT")         # indel
    assert not builder._usable_marker_alleles("A", "T")          # A/T palindrome
    assert not builder._usable_marker_alleles("C", "G")          # C/G palindrome
    assert not builder._usable_marker_alleles("A", "A")          # not a variant


def test_first_present_is_case_insensitive() -> None:
    row = {"Sample": "NA12345", "POP": "GBR"}
    assert builder._first_present(row, "sample_id", "sample") == "NA12345"
    assert builder._first_present(row, "Population", "pop") == "GBR"
    assert builder._first_present(row, "missing", "absent") == ""


def test_write_compact_panel_produces_expected_files(tmp_path: Path) -> None:
    samples = [
        {"sample_id": "S1", "population": "GBR", "superpopulation": "EUR", "sex": "male"},
        {"sample_id": "S2", "population": "YRI", "superpopulation": "AFR", "sex": "female"},
        {"sample_id": "S3", "population": "CHB", "superpopulation": "EAS", "sex": "female"},
        {"sample_id": "S4", "population": "PUR", "superpopulation": "AMR", "sex": "male"},
        {"sample_id": "S5", "population": "GIH", "superpopulation": "SAS", "sex": "female"},
    ]
    markers = [
        {"chrom": "1", "pos": 1000, "ref": "A", "alt": "G"},
        {"chrom": "2", "pos": 5000, "ref": "C", "alt": "T"},
        {"chrom": "3", "pos": 9000, "ref": "G", "alt": "A"},
    ]
    genotype_rows = [
        [0.0, 1.0, 2.0, 1.0, 0.0],
        [2.0, 0.0, 1.0, 1.0, 2.0],
        [1.0, 1.0, 0.0, 2.0, 1.0],
    ]
    result = builder.write_compact_panel(
        tmp_path,
        samples=samples,
        markers=markers,
        genotype_rows=genotype_rows,
        component_count=3,
        force=True,
    )
    assert result["status"] == "completed"
    for name in naming.PANEL_FILES:
        assert (tmp_path / name).exists(), name
    manifest = json.loads((tmp_path / naming.MANIFEST_NAME).read_text())
    assert manifest["panel_id"] == naming.PANEL_ID
    assert manifest["sample_count"] == 5
    assert manifest["marker_count"] == 3
    assert manifest["component_count"] == 3


def test_write_compact_panel_validates_shapes(tmp_path: Path) -> None:
    samples = [{"sample_id": "S1"}, {"sample_id": "S2"}]
    markers = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]
    with pytest.raises(ValueError):
        builder.write_compact_panel(
            tmp_path,
            samples=samples,
            markers=markers,
            genotype_rows=[[1.0]],         # wrong number of sample dosages
            force=True,
        )


def test_copy_panel_round_trips(tmp_path: Path) -> None:
    src = tmp_path / "src"
    dst = tmp_path / "dst"
    builder.write_compact_panel(
        src,
        samples=[{"sample_id": "S1", "population": "GBR", "superpopulation": "EUR", "sex": ""},
                 {"sample_id": "S2", "population": "YRI", "superpopulation": "AFR", "sex": ""},
                 {"sample_id": "S3", "population": "CHB", "superpopulation": "EAS", "sex": ""}],
        markers=[{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"},
                 {"chrom": "2", "pos": 200, "ref": "C", "alt": "T"}],
        genotype_rows=[[0.0, 1.0, 2.0], [2.0, 1.0, 0.0]],
        component_count=2,
        force=True,
    )
    builder.copy_panel(src, dst, force=True)
    for name in naming.PANEL_FILES:
        assert (dst / name).read_bytes() == (src / name).read_bytes()
