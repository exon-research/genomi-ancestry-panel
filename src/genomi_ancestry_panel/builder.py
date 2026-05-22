"""Build the compact 1000 Genomes 30x GRCh38 PCA reference panel.

The build streams all 22 autosomal phased VCFs from EBI, filters down to a
target marker count under MAF / call-rate / spacing constraints, decomposes
the resulting 3,202 × ~20,000 dosage matrix with SVD, and writes a TSV/JSON
panel directory ready for downstream PCA projection. See docs/math.md and
docs/filters.md for the underlying reasoning.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import math
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator

import numpy as np

from . import naming, source_urls


JsonObject = dict[str, Any]

DEFAULT_COMPONENT_COUNT = 10
DEFAULT_TARGET_MARKERS = 20_000
DEFAULT_MIN_SPACING_BP = 250_000
DEFAULT_MIN_CALL_RATE = 0.98
DEFAULT_MIN_MAF = 0.05
DEFAULT_MAX_MAF = 0.45

AUTOSOMES = source_urls.AUTOSOMES
BASES = {"A", "C", "G", "T"}


def build_1000g_30x_grch38_panel(
    *,
    output_dir: str | Path,
    force: bool = False,
    source_vcfs: Iterable[str | Path] | None = None,
    sample_metadata: str | Path | None = None,
    target_markers: int = DEFAULT_TARGET_MARKERS,
    min_spacing_bp: int = DEFAULT_MIN_SPACING_BP,
    component_count: int = DEFAULT_COMPONENT_COUNT,
) -> JsonObject:
    """End-to-end build: stream public VCFs, select markers, run PCA, write panel."""
    output = Path(output_dir)
    manifest = output / naming.MANIFEST_NAME
    if manifest.exists() and not force:
        return {"status": "cached", "library": naming.PANEL_LIBRARY, "manifest_path": str(manifest)}
    vcf_inputs = list(source_vcfs or source_urls.default_phased_vcf_urls())
    metadata_source = str(sample_metadata or source_urls.SAMPLE_PANEL_URL)
    samples_by_id = _load_sample_metadata(metadata_source)
    selected = _select_marker_genotypes(
        vcf_inputs,
        samples_by_id=samples_by_id,
        target_markers=target_markers,
        min_spacing_bp=min_spacing_bp,
    )
    if not selected["markers"]:
        raise RuntimeError("No ancestry panel markers passed deterministic filters.")
    return write_compact_panel(
        output,
        samples=selected["samples"],
        markers=selected["markers"],
        genotype_rows=selected["genotype_rows"],
        source_inputs={"sample_metadata": metadata_source, "vcfs": [str(value) for value in vcf_inputs]},
        component_count=component_count,
        force=True,
    )


def write_compact_panel(
    output_dir: str | Path,
    *,
    samples: list[JsonObject],
    markers: list[JsonObject],
    genotype_rows: list[list[float | None]],
    source_inputs: JsonObject | None = None,
    component_count: int = DEFAULT_COMPONENT_COUNT,
    force: bool = False,
) -> JsonObject:
    output = Path(output_dir)
    manifest_path = output / naming.MANIFEST_NAME
    if manifest_path.exists() and not force:
        return {"status": "cached", "library": naming.PANEL_LIBRARY, "manifest_path": str(manifest_path)}
    output.mkdir(parents=True, exist_ok=True)
    matrix = np.asarray(
        [[np.nan if value is None else float(value) for value in row] for row in genotype_rows],
        dtype=float,
    )
    if matrix.shape[0] != len(markers):
        raise ValueError("genotype_rows must have one row per marker")
    if matrix.shape[1] != len(samples):
        raise ValueError("genotype_rows must have one value per sample")
    means = np.nanmean(matrix, axis=1)
    scales = np.nanstd(matrix, axis=1, ddof=0)
    scales[~np.isfinite(scales) | (scales <= 0)] = 1.0
    imputed = np.where(np.isnan(matrix), means[:, None], matrix)
    standardized = ((imputed - means[:, None]) / scales[:, None]).T
    max_components = max(1, min(int(component_count), standardized.shape[0], standardized.shape[1]))
    u, singular_values, vt = np.linalg.svd(standardized, full_matrices=False)
    scores = u[:, :max_components] * singular_values[:max_components]
    loadings = vt[:max_components, :].T
    component_names = [f"PC{index + 1}" for index in range(max_components)]

    marker_rows = []
    for index, marker in enumerate(markers):
        marker_rows.append(
            {
                **marker,
                "marker_id": marker.get("marker_id") or f"{marker['chrom']}:{marker['pos']}:{marker['ref']}:{marker['alt']}",
                "mean": f"{float(means[index]):.10g}",
                "scale": f"{float(scales[index]):.10g}",
            }
        )

    _write_tsv(output / naming.SAMPLES_NAME, ["sample_id", "population", "superpopulation", "sex"], samples)
    _write_tsv(output / naming.MARKERS_NAME, ["marker_id", "chrom", "pos", "ref", "alt", "mean", "scale"], marker_rows)
    _write_tsv(
        output / naming.LOADINGS_NAME,
        ["marker_id", *component_names],
        [
            {
                "marker_id": marker_rows[index]["marker_id"],
                **{name: f"{float(loadings[index, pc_index]):.10g}" for pc_index, name in enumerate(component_names)},
            }
            for index in range(len(marker_rows))
        ],
    )
    _write_tsv(
        output / naming.REFERENCE_SCORES_NAME,
        ["sample_id", "population", "superpopulation", *component_names],
        [
            {
                "sample_id": sample["sample_id"],
                "population": sample.get("population", ""),
                "superpopulation": sample.get("superpopulation", ""),
                **{name: f"{float(scores[sample_index, pc_index]):.10g}" for pc_index, name in enumerate(component_names)},
            }
            for sample_index, sample in enumerate(samples)
        ],
    )
    stats = {
        "schema": "genomi-ancestry-panel-stats-v1",
        "sample_count": len(samples),
        "marker_count": len(marker_rows),
        "component_count": len(component_names),
        "target_marker_count": len(marker_rows),
        "built_at": _now(),
    }
    (output / naming.PANEL_STATS_NAME).write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "schema": "genomi-ancestry-reference-panel-v1",
        "panel_id": naming.PANEL_ID,
        "title": naming.PANEL_TITLE,
        "library": naming.PANEL_LIBRARY,
        "genome_build": "GRCh38",
        "method": "PCA projection from mean/std-normalized autosomal biallelic SNP dosages",
        "sample_count": len(samples),
        "marker_count": len(marker_rows),
        "component_count": len(component_names),
        "files": {
            "samples": naming.SAMPLES_NAME,
            "markers": naming.MARKERS_NAME,
            "pca_loadings": naming.LOADINGS_NAME,
            "reference_scores": naming.REFERENCE_SCORES_NAME,
            "panel_stats": naming.PANEL_STATS_NAME,
        },
        "source_urls": source_urls.source_urls(),
        "source_inputs": source_inputs or {},
        "built_at": _now(),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "status": "completed",
        "library": naming.PANEL_LIBRARY,
        "manifest_path": str(manifest_path),
        "panel_dir": str(output),
        "sample_count": len(samples),
        "marker_count": len(marker_rows),
        "component_count": len(component_names),
    }


def _select_marker_genotypes(
    sources: list[str | Path],
    *,
    samples_by_id: dict[str, JsonObject],
    target_markers: int,
    min_spacing_bp: int,
) -> JsonObject:
    selected_markers: list[JsonObject] = []
    genotype_rows: list[list[float | None]] = []
    samples: list[JsonObject] | None = None
    per_chrom_target = max(1, math.ceil(target_markers / len(AUTOSOMES)))
    selected_by_chrom = {chrom: 0 for chrom in AUTOSOMES}
    last_pos_by_chrom = {chrom: -min_spacing_bp for chrom in AUTOSOMES}
    for source in sources:
        header_samples: list[str] | None = None
        for line in _iter_text_lines(source):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_samples = line.rstrip("\n").split("\t")[9:]
                if samples is None:
                    samples = [_sample_record(sample_id, samples_by_id.get(sample_id, {})) for sample_id in header_samples]
                continue
            if header_samples is None or samples is None:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom = fields[0].removeprefix("chr")
            if chrom not in selected_by_chrom or selected_by_chrom[chrom] >= per_chrom_target:
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            if pos - last_pos_by_chrom[chrom] < min_spacing_bp:
                continue
            ref = fields[3].upper()
            alt = fields[4].split(",", 1)[0].upper()
            if not _usable_marker_alleles(ref, alt):
                continue
            dosages = _dosages_from_samples(fields[8], fields[9:])
            called = [value for value in dosages if value is not None]
            if not called or len(called) / len(dosages) < DEFAULT_MIN_CALL_RATE:
                continue
            alt_frequency = sum(called) / (2 * len(called))
            maf = min(alt_frequency, 1 - alt_frequency)
            if maf < DEFAULT_MIN_MAF or maf > DEFAULT_MAX_MAF:
                continue
            marker_id = fields[2] if fields[2] and fields[2] != "." else f"{chrom}:{pos}:{ref}:{alt}"
            selected_markers.append({"marker_id": marker_id, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt})
            genotype_rows.append(dosages)
            selected_by_chrom[chrom] += 1
            last_pos_by_chrom[chrom] = pos
            if len(selected_markers) >= target_markers:
                break
        if len(selected_markers) >= target_markers:
            break
    return {"samples": samples or [], "markers": selected_markers, "genotype_rows": genotype_rows}


def _load_sample_metadata(source: str | Path) -> dict[str, JsonObject]:
    try:
        lines = list(_iter_text_lines(source))
    except OSError:
        return {}
    if not lines:
        return {}
    reader = csv.DictReader(io.StringIO("".join(lines)), delimiter="\t")
    output: dict[str, JsonObject] = {}
    for row in reader:
        sample_id = _first_present(row, "SampleID", "sample", "Sample", "sample_id", "Individual ID")
        if not sample_id:
            continue
        output[sample_id] = {
            "population": _first_present(row, "Population", "population", "Population code", "pop"),
            "superpopulation": _first_present(row, "Superpopulation", "Superpopulation code", "superpopulation", "super_pop", "superpop"),
            "sex": _first_present(row, "Sex", "sex", "gender"),
        }
    return output


def _sample_record(sample_id: str, metadata: JsonObject) -> JsonObject:
    return {
        "sample_id": sample_id,
        "population": metadata.get("population") or "unknown",
        "superpopulation": metadata.get("superpopulation") or "unknown",
        "sex": metadata.get("sex") or "",
    }


def _first_present(row: dict[str, Any], *keys: str) -> str:
    lower = {key.lower(): value for key, value in row.items()}
    for key in keys:
        value = row.get(key)
        if value not in (None, ""):
            return str(value)
        value = lower.get(key.lower())
        if value not in (None, ""):
            return str(value)
    return ""


def _dosages_from_samples(format_field: str, sample_fields: list[str]) -> list[float | None]:
    keys = format_field.split(":")
    try:
        gt_index = keys.index("GT")
    except ValueError:
        return [None for _ in sample_fields]
    output: list[float | None] = []
    for sample in sample_fields:
        values = sample.split(":")
        if gt_index >= len(values):
            output.append(None)
            continue
        genotype = values[gt_index]
        if not genotype or "." in genotype:
            output.append(None)
            continue
        tokens = genotype.replace("|", "/").split("/")
        if any(token not in {"0", "1"} for token in tokens):
            output.append(None)
            continue
        output.append(float(sum(1 for token in tokens if token == "1")))
    return output


def _usable_marker_alleles(ref: str, alt: str) -> bool:
    if len(ref) != 1 or len(alt) != 1:
        return False
    if ref not in BASES or alt not in BASES or ref == alt:
        return False
    return {ref, alt} not in ({"A", "T"}, {"C", "G"})


def _iter_text_lines(source: str | Path) -> Iterator[str]:
    text_source = str(source)
    if text_source.startswith(("http://", "https://")):
        request = urllib.request.Request(text_source, headers={"User-Agent": "genomi-ancestry-panel/1.0"})
        with urllib.request.urlopen(request, timeout=120) as response:
            if text_source.endswith(".gz"):
                with gzip.GzipFile(fileobj=response) as gz_handle:
                    for raw in gz_handle:
                        yield raw.decode("utf-8", errors="replace")
            else:
                for raw in response:
                    yield raw.decode("utf-8", errors="replace")
        return
    path = Path(text_source).expanduser()
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        yield from handle


def _write_tsv(path: Path, fieldnames: list[str], rows: Iterable[JsonObject]) -> None:
    tmp = path.with_name(path.name + ".partial")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    tmp.replace(path)


def copy_panel(source_dir: str | Path, output_dir: str | Path, *, force: bool = False) -> JsonObject:
    """Copy an already-built panel directory into a target output dir."""
    source = Path(source_dir).expanduser()
    output = Path(output_dir).expanduser()
    if output.exists() and force:
        shutil.rmtree(output)
    output.mkdir(parents=True, exist_ok=True)
    for name in naming.PANEL_FILES:
        shutil.copyfile(source / name, output / name)
    return {"status": "completed", "library": naming.PANEL_LIBRARY, "manifest_path": str(output / naming.MANIFEST_NAME)}


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
