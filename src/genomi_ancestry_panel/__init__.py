"""Reproducible 1000 Genomes 30x GRCh38 PCA reference-panel builder.

This package is a self-contained pipeline that turns the public 1000 Genomes
30x release into a compact PCA reference panel for ancestry projection. The
output is a directory of TSV/JSON files (a few megabytes total) that
downstream tools can use to project a new sample's genotypes into the same
PC space as the 3,202 reference samples and look up nearest-neighbor
population labels.

Build with: `python -m genomi_ancestry_panel build --output panel-output/`
"""

from . import naming, source_urls
from .builder import (
    build_1000g_30x_grch38_panel,
    copy_panel,
    write_compact_panel,
    DEFAULT_COMPONENT_COUNT,
    DEFAULT_MAX_MAF,
    DEFAULT_MIN_CALL_RATE,
    DEFAULT_MIN_MAF,
    DEFAULT_MIN_SPACING_BP,
    DEFAULT_TARGET_MARKERS,
)

__all__ = [
    "build_1000g_30x_grch38_panel",
    "copy_panel",
    "write_compact_panel",
    "DEFAULT_COMPONENT_COUNT",
    "DEFAULT_MAX_MAF",
    "DEFAULT_MIN_CALL_RATE",
    "DEFAULT_MIN_MAF",
    "DEFAULT_MIN_SPACING_BP",
    "DEFAULT_TARGET_MARKERS",
    "naming",
    "source_urls",
]
