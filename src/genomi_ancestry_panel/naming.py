"""Stable panel identifiers and on-disk file names.

These names define the artifact contract. Downstream tools (Genomi and
anyone else consuming the panel) rely on them, so they should be considered
part of the public interface. Bump the panel version in
`pyproject.toml` if you change any name here.
"""

from __future__ import annotations

PANEL_ID = "1000g_30x_grch38"
PANEL_TITLE = "1000 Genomes 30x GRCh38 ancestry PCA panel"
PANEL_LIBRARY = "ancestry-1000g-30x-grch38"

MANIFEST_NAME = "manifest.json"
SAMPLES_NAME = "samples.tsv"
MARKERS_NAME = "markers.tsv"
LOADINGS_NAME = "pca_loadings.tsv"
REFERENCE_SCORES_NAME = "reference_scores.tsv"
PANEL_STATS_NAME = "panel_stats.json"

PANEL_FILES = (
    MANIFEST_NAME,
    SAMPLES_NAME,
    MARKERS_NAME,
    LOADINGS_NAME,
    REFERENCE_SCORES_NAME,
    PANEL_STATS_NAME,
)
