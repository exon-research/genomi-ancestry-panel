First public release of the compact 1000 Genomes 30x GRCh38 ancestry PCA
reference panel built by [`genomi-ancestry-panel`](https://github.com/exon-research/genomi-ancestry-panel).
The panel is a deterministic, ~3 MB summary of public 1000 Genomes data
that downstream tools (e.g. Genomi) use to project a new sample into a
shared principal-component space and look up nearest-neighbor reference
populations.

## Artifact

| File | Purpose |
| --- | --- |
| `panel-1000g-30x-grch38-1.0.0.tar.gz` | The full panel directory (samples, markers, loadings, reference scores, manifest, stats). |
| `SHA256SUMS` | SHA-256 checksum of the tarball. Verify before extracting. |

Verify and extract:

```bash
VERSION=1.0.0
curl -fsSL \
  "https://github.com/exon-research/genomi-ancestry-panel/releases/download/v${VERSION}/panel-1000g-30x-grch38-${VERSION}.tar.gz" \
  -o panel.tar.gz
sha256sum -c <(curl -fsSL "https://github.com/exon-research/genomi-ancestry-panel/releases/download/v${VERSION}/SHA256SUMS")
mkdir -p ~/.genomi/reference/ancestry
tar -xzf panel.tar.gz -C ~/.genomi/reference/ancestry
```

## What's inside

After extraction at `~/.genomi/reference/ancestry/1000g_30x_grch38/`:

- `manifest.json` — panel id, version, source URLs, build timestamp, file index
- `samples.tsv` — 3,202 reference samples with population / superpopulation / sex labels (2,504 phase-3 unrelated + 698 trio extensions; the latter carry `unknown` labels by design)
- `markers.tsv` — ~20,000 selected SNPs with `chrom / pos / ref / alt / per-marker mean / scale`
- `pca_loadings.tsv` — projection matrix (markers × principal components)
- `reference_scores.tsv` — coordinates of each reference sample in PC space
- `panel_stats.json` — counts and build metadata

The full 3,202 × ~20,000 dosage matrix that the builder computes in memory is **not** stored — only the PCA-derived summaries that are sufficient to project new samples.

## Source data

- **1000 Genomes 30x phased VCFs** at `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/`
- **Sample / population panel** at `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel`
- **Methods reference**: Byrska-Bishop et al. 2022 ([DOI](https://doi.org/10.1016/j.cell.2022.08.004))

All inputs are publicly available under the [1000 Genomes Project reuse policy](http://www.internationalgenome.org/IGSR_disclaimer).

## Build parameters

| Constant | Value | Why |
| --- | --- | --- |
| `DEFAULT_TARGET_MARKERS` | 20,000 | Comfortable margin over the ~5–10K independent markers needed for stable PC structure. |
| `DEFAULT_MIN_MAF` | 0.05 | Drop rare alleles that don't carry population-structure signal at N=3,202. |
| `DEFAULT_MAX_MAF` | 0.45 | Skip near-50/50 markers vulnerable to numerical edge cases. |
| `DEFAULT_MIN_CALL_RATE` | 0.98 | Keep markers where mean imputation doesn't measurably bias PC coordinates. |
| `DEFAULT_MIN_SPACING_BP` | 250,000 | Spacing-based LD prune so PC directions aren't double-weighted by haplotype-block redundancy. |
| `DEFAULT_COMPONENT_COUNT` | 10 | Captures continental + macro-regional structure; past 10 PCs is mostly noise. |

See [`docs/filters.md`](https://github.com/exon-research/genomi-ancestry-panel/blob/master/docs/filters.md) for justifications and citations.

## Reproducibility

This release was built from commit `<builder_commit_sha>` of the [`genomi-ancestry-panel`](https://github.com/exon-research/genomi-ancestry-panel) repo. Anyone can verify the artifact by re-running:

```bash
git clone https://github.com/exon-research/genomi-ancestry-panel.git
cd genomi-ancestry-panel && git checkout v1.0.0
uv venv .venv && uv pip install -e ".[dev]" --python .venv/bin/python
.venv/bin/python -m genomi_ancestry_panel build --output panel-output/
```

Expect ~3 hours of compute and ~300–500 GB of streamed (not stored) network traffic.

## What this panel is *not*

- **Not ethnicity, nationality, race, or identity prediction.** Reference-panel similarity in PC space is similarity to the 26 specific 1000 Genomes cohorts. It is not a determination of where you or your ancestors are from.
- **No admixture proportions.** This is single-point PCA projection, not unsupervised mixture-model deconvolution (ADMIXTURE / RFMix-style).
- **No haplogroups.** mtDNA and Y-chromosome lineages are inherited differently and require dedicated tools.
- **No local ancestry / chromosome painting.**
- **No relative matching / IBD.**
- **No archaic introgression** (Neanderthal / Denisovan).

These are deliberate scope limits, not bugs.

## License

MIT for this code and packaging. The underlying 1000 Genomes data retains its [public reuse policy](http://www.internationalgenome.org/IGSR_disclaimer). See [LICENSE](https://github.com/exon-research/genomi-ancestry-panel/blob/master/LICENSE) and [`docs/citations.md`](https://github.com/exon-research/genomi-ancestry-panel/blob/master/docs/citations.md) for full attribution.

If you use this artifact in published research, please cite:

1. Byrska-Bishop M et al. 2022 — for the 1000 Genomes 30x source data
2. Patterson, Price & Reich 2006 — for the PCA-for-ancestry methodology
3. This release tag: `exon-research/genomi-ancestry-panel@v1.0.0`
