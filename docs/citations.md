# Citations and source data

## 1000 Genomes Project

The reference samples and phased VCFs come from the 1000 Genomes Project's 30x re-sequencing collection.

- **30x re-release publication**: Byrska-Bishop M et al. 2022. *High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios.* Cell 185 (18): 3426–3440.e19. <https://doi.org/10.1016/j.cell.2022.08.004>
- **IGSR portal**: <https://www.internationalgenome.org/data-portal/data-collections/30x-grch38.html>
- **Working directory used by this builder**: <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/>
- **Phased VCFs used by this builder**: <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/>
- **Phase-3 sample / population panel used for sample labels**: <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel>
- **Original 1000 Genomes phase-3 paper**: The 1000 Genomes Project Consortium. 2015. *A global reference for human genetic variation.* Nature 526: 68–74. <https://doi.org/10.1038/nature15393>

All 1000 Genomes data are available under the project's public reuse policy: <http://www.internationalgenome.org/IGSR_disclaimer>.

## PCA for population structure

- Patterson N, Price AL, Reich D. 2006. *Population structure and eigenanalysis.* PLoS Genet 2(12):e190. <https://doi.org/10.1371/journal.pgen.0020190> — the canonical PCA-for-ancestry paper, including the per-marker normalization scheme this builder follows.
- Price AL et al. 2006. *Principal components analysis corrects for stratification in genome-wide association studies.* Nat Genet 38: 904–909. <https://doi.org/10.1038/ng1847>
- Conomos MP, Reiner AP, Weir BS, Thornton TA. 2015. *Robust inference of population structure for ancestry prediction and correction of stratification in the presence of relatedness.* Genet Epidemiol. <https://doi.org/10.1002/gepi.21896>
- Galinsky KJ et al. 2016. *Fast Principal-Component Analysis Reveals Convergent Evolution of ADH1B in Europe and East Asia.* AJHG 98:456–472. <https://doi.org/10.1016/j.ajhg.2015.12.022>

## Linkage disequilibrium and marker pruning

- Slatkin M. 2008. *Linkage disequilibrium — understanding the evolutionary past and mapping the medical future.* Nat Rev Genet 9: 477–485.
- The 1000 Genomes Project Consortium. 2010. *A map of human genome variation from population-scale sequencing.* Nature 467: 1061–1073.

## Tooling that informed this design

- **PLINK 2**: Chang CC et al. 2015. <https://doi.org/10.1186/s13742-015-0047-8> — `--indep-pairwise` LD pruning and `--pca` for the standard reference implementation.
- **SNPRelate / GENESIS**: Zheng X et al. 2012. <https://doi.org/10.1093/bioinformatics/bts606> — efficient R implementation of population PCA.

## How to cite this artifact

If you use a panel artifact built from this repo in published research, please cite:

1. Byrska-Bishop et al. 2022 (above) for the 1000 Genomes 30x data.
2. Patterson, Price & Reich 2006 (above) for the PCA-for-ancestry methodology.
3. This repository at the published release tag (`genomi-ancestry-panel@vX.Y.Z`).
