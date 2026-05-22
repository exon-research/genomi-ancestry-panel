# genomi-ancestry-panel

A reproducible builder for a **compact 1000 Genomes 30x GRCh38 PCA reference panel**, plus the math and reasoning behind it.

The output is a ~3 MB directory of TSV/JSON files that downstream tools (e.g. Genomi) consume to project a new sample's genotypes into the same PC space as 3,202 public 1000 Genomes reference samples, and look up nearest-neighbor population labels.

## How to use the built panel

You don't normally interact with this repo directly. Downstream tools download the tarball from the matching GitHub release, extract it under their reference-data root, and use it for PCA projection.

Manual download (for any consumer):

```bash
VERSION=1.0.0
mkdir -p ~/.genomi/reference/ancestry/
curl -fsSL \
  "https://github.com/exon-research/genomi-ancestry-panel/releases/download/v${VERSION}/panel-1000g-30x-grch38-${VERSION}.tar.gz" \
  | tar -xz -C ~/.genomi/reference/ancestry/
```

The extracted directory is what Genomi's `ancestry.estimate_population_context` reads at query time. Projection is seconds, not hours.

## How to rebuild the panel from scratch

You should only need this if you are publishing a new release, or auditing the build.

```bash
git clone https://github.com/exon-research/genomi-ancestry-panel.git
cd genomi-ancestry-panel
uv venv .venv && uv pip install -e ".[dev]" --python .venv/bin/python
.venv/bin/python -m genomi_ancestry_panel build --output panel-output/
```

Expect:

- **Time:** ~3 hours on a typical VM (single-threaded Python, network-bound by EBI).
- **Network:** streams roughly 300–500 GB of compressed VCF data (does not save to disk).
- **Memory:** peaks at ~400 MB.
- **Disk:** ~3 MB final.

The build is deterministic given the same upstream sources and the same constants in [`builder.py`](src/genomi_ancestry_panel/builder.py).

## How to cut a release

```bash
# 1. bump version in pyproject.toml
# 2. build + tarball + sha256:
./scripts/build_and_release.sh
# 3. attach to GitHub release:
gh release create vX.Y.Z \
  release-staging/panel-1000g-30x-grch38-X.Y.Z.tar.gz \
  release-staging/SHA256SUMS \
  --title "Panel vX.Y.Z" \
  --notes-file CHANGELOG.md
```

## The math and the reasoning, in detail

- [docs/pca_intuition.md](docs/pca_intuition.md) — Plain-English walkthrough: SNPs, genotypes, MAF, why PCA picks out continental structure.
- [docs/math.md](docs/math.md) — Formal: SVD-based PCA derivation, normalization, why projection of a new sample reduces to a vector–matrix product.
- [docs/filters.md](docs/filters.md) — Justifications for `DEFAULT_MIN_MAF=0.05`, `DEFAULT_MAX_MAF=0.45`, `DEFAULT_MIN_CALL_RATE=0.98`, `DEFAULT_MIN_SPACING_BP=250000`, `DEFAULT_TARGET_MARKERS=20000`, `DEFAULT_COMPONENT_COUNT=10`. With citations.
- [docs/citations.md](docs/citations.md) — 1000 Genomes Project references, McVean & Pritchard methods papers, PCA-for-ancestry literature.

## What's in the panel artifact

A single tarball, e.g. `panel-1000g-30x-grch38-1.0.0.tar.gz`, that extracts to:

```
panel-1000g-30x-grch38-1.0.0/
├── manifest.json           # panel id, version, source URLs, build timestamp, file index
├── samples.tsv             # 3,202 reference samples with population / superpopulation / sex labels
├── markers.tsv             # ~20,000 selected SNPs with chrom / pos / ref / alt / per-marker mean+scale
├── pca_loadings.tsv        # the projection matrix (markers × principal components)
├── reference_scores.tsv    # where each reference sample sits in PC space
└── panel_stats.json        # counts, build metadata
```

Total disk footprint: a few megabytes. The full 3,202 × ~20,000 dosage matrix that the build computes in memory is *not* stored — only the PCA-derived summaries that are sufficient to project new samples.

## What the panel does *not* do

- **Not ethnicity, nationality, race, or identity prediction.** Reference-panel similarity in PC space is exactly that: similarity to the 26 specific 1000 Genomes cohorts. It is not a determination of where you or your ancestors are from.
- **No admixture proportions.** This is single-point PCA projection, not unsupervised mixture-model deconvolution (ADMIXTURE / RFMix-style).
- **No haplogroups.** mtDNA and Y-chromosome lineages are inherited differently and require dedicated tools.
- **No local ancestry / chromosome painting.** This panel returns a single point in PC space per sample, not per-locus ancestry calls.
- **No relative matching / IBD.** No relative search.
- **No archaic introgression.** Neanderthal / Denisovan ancestry estimation needs a separate archaic-reference dataset.

These are deliberate scope limits, not bugs.

## Provenance and reproducibility

Every panel artifact published from this repo carries a `manifest.json` that records:

- `panel_id`, `library`, `genome_build`
- `source_urls` — exact public 1000 Genomes URLs streamed
- `source_inputs.vcfs` — the 22 chromosome VCF URLs (or local overrides)
- `sample_count`, `marker_count`, `component_count`
- `built_at` — UTC timestamp
- (recommended for future releases) `builder_commit_sha` and per-file SHA256s

Anyone can rebuild from this repo's source and confirm they get the same output, modulo nondeterministic dict iteration that doesn't affect the final loadings/scores.

## License

MIT, except for the 1000 Genomes data which retains its [public reuse policy](http://www.internationalgenome.org/IGSR_disclaimer). See [LICENSE](LICENSE) and [docs/citations.md](docs/citations.md).
