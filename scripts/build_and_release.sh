#!/usr/bin/env bash
# Build the panel, tarball it, generate SHA256SUMS, leave artifacts in release-staging/.
# Upload the tarball + SHA256SUMS to the matching GitHub release tag manually:
#     gh release create v$VERSION release-staging/panel-1000g-30x-grch38-$VERSION.tar.gz release-staging/SHA256SUMS
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VERSION="$(python3 -c "import tomllib,pathlib; print(tomllib.loads(pathlib.Path('$ROOT/pyproject.toml').read_text())['project']['version'])")"
STAGE="$ROOT/release-staging"
OUT="$STAGE/panel-1000g-30x-grch38-$VERSION"
TARBALL="$STAGE/panel-1000g-30x-grch38-$VERSION.tar.gz"

rm -rf "$STAGE"
mkdir -p "$STAGE"

echo "==> Building panel into $OUT (this takes ~3 hours and streams hundreds of GB)"
python3 -m genomi_ancestry_panel build --output "$OUT"

echo "==> Creating tarball $TARBALL"
tar -czf "$TARBALL" -C "$STAGE" "$(basename "$OUT")"

echo "==> Computing SHA256SUMS"
( cd "$STAGE" && sha256sum "$(basename "$TARBALL")" > SHA256SUMS )

echo "==> Tarball ready:"
ls -lh "$TARBALL" "$STAGE/SHA256SUMS"
echo ""
echo "Next steps:"
echo "  gh release create v$VERSION $TARBALL $STAGE/SHA256SUMS --title 'Panel v$VERSION'"
