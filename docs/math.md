# The math, formally

## Inputs

- A set of $N = 3{,}202$ phased reference samples from the 1000 Genomes 30x release.
- A set of $M$ ancestry-informative biallelic SNPs (target $M \approx 20{,}000$ after filtering).
- A dosage matrix $X \in \{0, 1, 2, \text{NaN}\}^{M \times N}$ where $X_{m,n}$ is the count of ALT alleles for sample $n$ at marker $m$ (or NaN if not called).

## Per-marker normalization

For each marker $m$:

- Compute mean dosage $\mu_m = \operatorname{nanmean}(X_{m,\cdot})$.
- Compute standard deviation $\sigma_m = \operatorname{nanstd}(X_{m,\cdot}, \text{ddof}=0)$.
- If $\sigma_m \le 0$ or is non-finite, set $\sigma_m = 1$ (treat as degenerate-but-keep).
- Replace NaNs in row $m$ with $\mu_m$ (mean imputation).

Define the standardized matrix:

$$Z_{m,n} = \frac{X_{m,n} - \mu_m}{\sigma_m}.$$

Then transpose so rows are samples and columns are markers:

$$A = Z^\top \in \mathbb{R}^{N \times M}.$$

This per-marker centering is the standard approach for population PCA (see Patterson, Price & Reich 2006). It removes the global allele-frequency level from each marker so PCA picks up *differences* in frequency between population structure rather than overall mean dosage.

## SVD

Compute the thin SVD:

$$A = U \, S \, V^\top, \quad U \in \mathbb{R}^{N \times k}, \; S \in \mathbb{R}^{k \times k}, \; V \in \mathbb{R}^{M \times k},$$

where $k = \min(N, M)$.

Choose the number of principal components to retain, $K$ (default 10). Then:

- **Reference scores** in PC space: $\Phi = U_{:, 1{:}K} \cdot \operatorname{diag}(S_{1{:}K})$, with shape $N \times K$. Row $n$ is sample $n$'s coordinates.
- **Loadings**: $L = V_{:, 1{:}K}$, with shape $M \times K$. Row $m$ is marker $m$'s contribution to each PC.

These are the only two large arrays written to disk. The original dosage matrix is discarded.

## Projecting a new sample

Given a new sample with genotype dosages $x \in \mathbb{R}^M$ at the panel marker positions (with missing values imputed using the panel's $\mu_m$):

1. Standardize: $z_m = (x_m - \mu_m) / \sigma_m$.
2. Project: $\phi = z^\top L$, with shape $K$.
3. The new sample's coordinates in PC space are $\phi$.

This is a single vector–matrix product. It does *not* require redoing SVD, and it does not require reading the full reference dosage matrix. That's the whole point of caching the panel as loadings + reference scores instead of as raw genotypes.

## Nearest-neighbor lookup

Once $\phi$ is computed:

- For each reference sample $n$, compute Euclidean distance $d_n = \lVert \phi - \Phi_{n,\cdot} \rVert_2$.
- Sort by distance.
- Report the top-$K_{nn}$ nearest reference samples and their population labels.
- Optionally aggregate by population/superpopulation: the **nearest population centroid** is the centroid in PC space that minimizes Euclidean distance to $\phi$.

The downstream Genomi tool (`ancestry.estimate_population_context`) does this lookup.

## Why this stops being meaningful past ~30 PCs

Population structure in modern humans is concentrated in a small number of dimensions because:

- The deepest split (Africa vs. non-Africa) is a single direction.
- The major continental splits are a handful more.
- Below that, regional subpopulation differences add a few PCs.
- Past ~30 PCs you start fitting individual-level noise (recent ancestry, sequencing batch effects, IBD chunks) that doesn't generalize.

Empirically, retaining 10 PCs captures essentially all the major continental and macro-regional signal. Going higher gives diminishing returns and risks overfitting. (See Patterson et al. 2006; Galinsky et al. 2016 for fastPCA convergence discussion.)

## Determinism notes

- The SVD itself is deterministic given the same input matrix.
- Sign of each PC is ambiguous (you can flip the sign of column $L_{:, k}$ and column $\Phi_{:, k}$ simultaneously without changing the geometry); some BLAS implementations flip signs nondeterministically. Output coordinates may differ in sign across runs, but distances and nearest-neighbor structure are preserved.
- Marker selection iterates VCF lines in deterministic order; the only nondeterminism is whatever the underlying SVD library chooses, which doesn't affect downstream nearest-neighbor lookup.
