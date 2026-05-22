# PCA-for-ancestry, in plain English

A guide for clinicians, journalists, and curious users. No prior genetics or linear algebra assumed.

## DNA, in one paragraph

Your DNA is a chain of ~3 billion letters from {A, C, G, T} (called **bases**), broken into 23 pieces called **chromosomes**. You have two copies of each non-sex chromosome, one from each parent. At ~99.9% of positions every human is identical. At ~100 million positions some humans carry a different letter than others — those positions are catalogued as **SNPs** ("snips," for Single Nucleotide Polymorphisms). A SNP has a fixed position and two letters that exist at it: **REF** (in the reference genome) and **ALT** (the variant). Your **genotype** at a SNP is the pair of letters you specifically have — REF/REF, REF/ALT, or ALT/ALT, encoded `0/0`, `0/1`, `1/1`.

## Why some SNPs are useful for ancestry

When humans spread out of Africa over the last ~70,000 years, populations got geographically isolated. Random genetic drift and local selection nudged allele frequencies in different directions in different places. So today, some SNPs have very different ALT frequencies across populations — e.g. an allele common in northern Europeans is rare in West Africans. Aggregating signal across many such SNPs lets us reconstruct broad continental groupings.

The SNPs that carry the most ancestry signal are the **common ones** — minor allele frequency (MAF) in the 5%–50% range. Rare alleles are usually private to families and tell you nothing useful about populations; very common ones in *everyone* don't separate groups. So you filter to the informative middle.

## What PCA does

Imagine your sample as a row of 20,000 numbers — your dosages (0, 1, or 2) at each ancestry-informative SNP. Imagine each of the 3,202 reference samples as the same row. You now have a 3,202 × 20,000 matrix.

**PCA** ("Principal Component Analysis") finds the directions of greatest variation in that matrix. PC1 is the single axis along which the 3,202 reference samples spread out the most. PC2 is the next-best independent direction. PC3, PC4, and so on capture finer structure.

For human variation it turns out:

- **PC1** separates Africans from non-Africans — the deepest split in human ancestry.
- **PC2** separates East Asians from Europeans and South Asians.
- **PC3** separates South Asians from Europeans.
- **PC4–10** capture regional structure within continents.

So a sample's full ancestry signal compresses down to ~10 numbers (its coordinates on the first 10 PCs). This is the compression that makes the panel small.

## Projecting a new sample

Once the panel exists, projecting *your* genotypes is fast:

1. Read your genotype at each of the 20,000 panel marker positions from your VCF.
2. Convert each genotype to a dosage (`0/0 → 0`, `0/1 → 1`, `1/1 → 2`). Missing positions get the panel's per-marker mean.
3. Center and scale by the panel's per-marker mean and standard deviation.
4. Multiply by the loadings matrix (`markers × PCs`) → your coordinates in PC space.
5. Compute distance to all 3,202 reference samples in PC space.
6. Aggregate the nearest neighbors' population labels → "you sit closest to GBR / CEU / IBS / FIN / PUR samples."

Steps 1–5 are seconds of CPU. The slow work was the panel build, and that's done once for everyone.

## What this gives you (and what it doesn't)

The output is "your sample sits in this region of a low-dim space defined by 3,202 reference samples from 26 specific populations." That's similarity, not identity.

It is **not**:

- An ethnicity label.
- A nationality.
- A determination of where your ancestors came from.
- An admixture pie chart (Genomi panel does not deconvolve mixed ancestry).
- A relative-matching result.
- Coverage of every world population (Indigenous Americas, Pacific Islanders, Middle East, and large parts of Africa and Central Asia are underrepresented in 1000 Genomes).

It is **really good for**:

- "Does this sample broadly look European / East Asian / South Asian / African / Admixed American?"
- "Which 1000G reference cohort sits closest in PC space?"
- Sanity-checking that a sequencing run is what you think it is (e.g. flagging mislabelled samples).

For everything beyond that (admixture %, family lineage, archaic ancestry), use dedicated tools — the limits are baked into the method, not a Genomi shortcoming.
