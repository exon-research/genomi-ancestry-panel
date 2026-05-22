"""Command-line entry point for the ancestry panel builder.

    python -m genomi_ancestry_panel build --output panel-output/
    python -m genomi_ancestry_panel copy --from /existing/panel --to /target/panel
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import builder


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="genomi-ancestry-panel", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    build_parser = sub.add_parser("build", help="Build the 1000G 30x GRCh38 PCA panel from public source data.")
    build_parser.add_argument("--output", required=True, help="Directory to write the panel files into.")
    build_parser.add_argument("--force", action="store_true", help="Rebuild even if a manifest already exists.")
    build_parser.add_argument("--target-markers", type=int, default=builder.DEFAULT_TARGET_MARKERS)
    build_parser.add_argument("--min-spacing-bp", type=int, default=builder.DEFAULT_MIN_SPACING_BP)
    build_parser.add_argument("--component-count", type=int, default=builder.DEFAULT_COMPONENT_COUNT)
    build_parser.add_argument("--source-vcfs", nargs="*", help="Override phased VCF sources (URL or local path).")
    build_parser.add_argument("--sample-metadata", help="Override sample/population panel file source.")

    copy_parser = sub.add_parser("copy", help="Copy an already-built panel directory into a new location.")
    copy_parser.add_argument("--from", dest="source", required=True)
    copy_parser.add_argument("--to", dest="target", required=True)
    copy_parser.add_argument("--force", action="store_true")

    args = parser.parse_args(argv)

    if args.command == "build":
        result = builder.build_1000g_30x_grch38_panel(
            output_dir=args.output,
            force=args.force,
            source_vcfs=args.source_vcfs,
            sample_metadata=args.sample_metadata,
            target_markers=args.target_markers,
            min_spacing_bp=args.min_spacing_bp,
            component_count=args.component_count,
        )
    elif args.command == "copy":
        result = builder.copy_panel(args.source, args.target, force=args.force)
    else:
        parser.error(f"Unknown command: {args.command}")
        return 2

    json.dump(result, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
