"""Public source data URLs the builder pulls from.

All inputs are public, license-clear, and stable per 1000 Genomes data
release. If the upstream paths change, update them here and bump the panel
version. The publication URL anchors the methodology to a citeable source.
"""

from __future__ import annotations

IGSR_COLLECTION_URL = "https://www.internationalgenome.org/data-portal/data-collections/30x-grch38.html"
IGSR_WORKING_DIR_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
PHASED_VCF_DIR_URL = IGSR_WORKING_DIR_URL + "20201028_3202_phased/"

# Phase-3 panel file carries sample / pop / super_pop / gender for the 2,504
# unrelated reference samples. The pedigree file in the working directory only
# carries sample / parent / sex, so it cannot label population on its own.
SAMPLE_PANEL_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

PUBLICATION_URL = "https://doi.org/10.1016/j.cell.2022.08.004"

AUTOSOMES = tuple(str(index) for index in range(1, 23))


def default_phased_vcf_urls() -> list[str]:
    return [
        f"{PHASED_VCF_DIR_URL}CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"
        for chrom in AUTOSOMES
    ]


def source_urls() -> dict[str, str]:
    return {
        "igsr_collection": IGSR_COLLECTION_URL,
        "working_directory": IGSR_WORKING_DIR_URL,
        "sample_panel": SAMPLE_PANEL_URL,
        "phased_vcf_directory": PHASED_VCF_DIR_URL,
        "publication": PUBLICATION_URL,
    }
