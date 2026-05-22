# Filter constants and their justifications

The builder filters the ~80 million phased variants in the 1000 Genomes 30x release down to ~20,000 ancestry-informative markers. Each filter is implemented in [`src/genomi_ancestry_panel/builder.py`](../src/genomi_ancestry_panel/builder.py). This document explains *why* each constant has its value.

Changing any of these constants is a versioned, reviewable decision — bump the panel version in `pyproject.toml` and document the change in a release note.

## `DEFAULT_MIN_MAF = 0.05` and `DEFAULT_MAX_MAF = 0.45`

Filter: keep markers with population minor allele frequency in `[0.05, 0.45]`.

**Why a lower bound at 0.05.** Below 5% MAF, an allele is rare enough that most of the 3,202 reference samples carry zero copies. PCA on a column with extreme skew (mostly zeros + a handful of 1s) over-weights individual outliers and doesn't reflect population structure. Practical lower bounds used in the literature range from 0.01 (when N is huge) to 0.05 (when N is in the few-thousand range). For 3,202 samples, 0.05 means roughly ≥300 ALT chromosomes — enough mass that the standardization step (subtracting $\mu_m$) is well-behaved.

**Why an upper bound at 0.45.** Strictly speaking MAF caps at 0.5 by definition. Capping at 0.45 instead is a small margin against numerical edge cases (markers right at 0.5 sometimes have weird normalization behavior). Has minimal effect on the marker pool because the distribution thins rapidly toward 0.5.

**Why this matters for ancestry.** Common variants (MAF 5%–45%) are old — they predate the out-of-Africa migration and have accumulated frequency differences across populations. Rare variants are mostly recent and population-private, useful for IBD/family-finding but not for PCA-based ancestry.

References:
- Patterson, Price, Reich. 2006. *Population structure and eigenanalysis.* PLoS Genet.
- Conomos et al. 2015. *Robust inference of population structure for ancestry prediction and correction of stratification in the presence of relatedness.* Genet Epidemiol.

## `DEFAULT_MIN_CALL_RATE = 0.98`

Filter: keep markers where at least 98% of the 3,202 reference samples have a confident genotype call.

**Why 0.98.** PCA cannot handle missing data directly — you have to impute or drop. The builder uses mean imputation (`X_{m,n} ← μ_m` when missing). Imputation introduces bias proportional to the missingness rate: replacing 5% of dosages with the column mean pulls those samples toward the column centroid, distorting PC coordinates. Capping missingness at 2% (call rate ≥ 0.98) keeps that bias small.

The 1000 Genomes 30x release is high-quality WGS; most variants in the phased set have call rates above 0.99 already, so this filter mainly excludes positions in difficult-to-call regions (centromeres, segmental duplications, low-mappability).

## `DEFAULT_MIN_SPACING_BP = 250000`

Filter: after selecting a marker at position `p` on chromosome `c`, ignore all candidate markers within 250 kb downstream on the same chromosome.

**Why spacing at all.** Adjacent SNPs are correlated — they're inherited together as chunks called **haplotype blocks**. This is **linkage disequilibrium (LD)**. Including 10 highly-correlated SNPs in PCA effectively gives that haplotype 10x the weight of an independent locus, distorting the PC basis. Cheap fix: keep only one marker per LD block.

**Why 250 kb specifically.** Human LD blocks are population-dependent but typically span 10–100 kb. Long-range LD can extend further in certain regions (the MHC on chr6, for example, can have LD blocks >1 Mb). A 250 kb spacing rule is a conservative middle ground that aggressively prunes most LD-correlated pairs without being so wide that we run out of candidate markers on smaller chromosomes.

Alternative implementations use explicit pairwise $r^2$ filtering (e.g. PLINK's `--indep-pairwise 50 5 0.2`). The spacing approximation is faster (no second pass), trades some precision for simplicity, and has the same first-order effect.

## `DEFAULT_TARGET_MARKERS = 20000`

Filter: stop selecting once we have 20,000 markers total (with a per-chromosome quota of `ceil(20000 / 22) ≈ 910`).

**Why 20,000.** Empirically, ancestry PCA on 3,202 samples stabilizes well below 10,000 informative markers. We pick 20,000 as a comfortable margin — high enough that variance in marker selection doesn't shift PC directions, low enough that build time stays bounded and the panel fits in single-digit MB.

Going higher gives diminishing returns because:
- Past ~10,000 markers, additional SNPs are correlated with existing ones (we've used up the largely-independent loci).
- PC directions are already well-resolved.

Going lower (under ~5,000) starts to introduce sampling noise into PC directions, particularly for the deeper PCs that capture sub-continental structure.

## `DEFAULT_COMPONENT_COUNT = 10`

Retain the first 10 principal components.

**Why 10.** For human population data:
- PC1: AFR vs. non-AFR
- PC2: EAS vs. EUR + SAS
- PC3: SAS vs. EUR
- PC4–7: regional structure within continents (e.g. N vs. S European, N vs. S East Asian)
- PC8–10: finer subpopulation structure and admixture axes

Beyond PC10, you start picking up batch-level noise rather than generalizable structure. Storing only 10 keeps the panel small (loadings file is `M × 10`, reference_scores file is `N × 10`).

A consumer that needs more dimensions can rebuild with `--component-count 30` and republish.

## Biallelic SNPs only

Filter: only single-letter REF and ALT, both in `{A, C, G, T}`, and not the strand-ambiguous pairs `{A, T}` or `{C, G}`.

**Why single-letter only.** Indels (insertions/deletions) are harder to genotype consistently and have lower call rates. They also don't have a clean dosage interpretation for PCA.

**Why drop A/T and C/G.** These are **strand-ambiguous palindromes**: if your VCF reports ALT=A and the reference panel reports ALT=T at the same SNP, you can't tell whether you and the panel agree (different strands of the same allele) or disagree (different alleles) without knowing the strand each used. Dropping them avoids a class of silent projection errors. Standard practice in population-genetics pipelines (PLINK, BCFtools, ADMIXTURE).

## Autosomes only

Filter: chromosomes 1–22. X, Y, and mitochondrial DNA are excluded.

**Why.** X is half-dose in males (one X copy) — different scale. Y and mitochondria are inherited from only one parent, follow different inheritance, and are typically used for haplogroup analysis rather than continental ancestry PCA. Including them would distort the math.

## Phase-3 sample panel for population labels

The builder loads sample-to-population mapping from `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel`. This is the 2,504-sample phase-3 panel file with `sample / pop / super_pop / gender` columns.

The 30x release includes ~698 additional trio extension samples (parents of phase-3 samples) that are not in the phase-3 panel file. Those samples appear in the VCF headers and get included in the PCA matrix, but their population label falls through to `unknown`.

This is the right behavior: the trio extensions are *related* to the unrelated phase-3 cohort and should not bias the PCA basis with their additional correlation. Their PC coordinates are computed, but they don't anchor any population centroid.
