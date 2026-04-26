# Intestinal, water, and sediment microbiota of Pacific white shrimp (*Litopenaeus vannamei*) ponds fed graded levels of corn fermented protein

16S rRNA amplicon sequencing analysis of microbial community structure across intestine, sediment, and water compartments in *Litopenaeus vannamei* production ponds fed graded dietary levels of corn fermented protein (CFP).

---

## Overview

This repository contains R scripts for reproducing the microbiome analyses reported in "Intestinal, water, and sediment microbiota of Pacific white shrimp (*Litopenaeus vannamei*) ponds fed graded levels of corn fermented protein." The study characterizes microbial diversity, composition, and differential abundance across three compartments — intestine, pond sediment, and pond water — in *L. vannamei* production ponds fed CFP at 0, 5, 10, and 20 g 100 g⁻¹ over a 10-week grow-out period.

Raw 16S rRNA amplicon sequencing data are deposited in the NCBI Sequence Read Archive under BioProject accession **PRJNA1445501**.

---

## Repository Structure

```
.
├── README.md
├── metadata.txt                  # Sample metadata (see Metadata section)
└── scripts/
    ├── 01_16S_Pipeline.R         # 16S processing, seqtab_nochim.rds → step4_phyloseq_object.rds
    ├── 02_Rarefaction.R          # Rarefaction curve analysis and depth selection
    ├── 03_ThresholdAnalysis.R    # Abundance and prevalence filtering threshold analysis
    ├── 04_AlphaDiversity.R       # Alpha diversity, richness, Shannon, Berger-Parker, Faith's PD
    ├── 05_BetaDiversity.R        # Beta diversity, PCoA, PERMANOVA, BETADISPER
    ├── 06_VennDiagram.R          # Genus sharing across dietary treatments within each compartment
    ├── 07_RelativeAbundance.R    # Relative abundance alluvial plots
    └── 08_DifferentialAnalysis.R # ANCOM-BC2 differential abundance analysis
```

---

## Data Availability

| Resource | Location |
|---|---|
| Raw FASTQ reads | NCBI SRA, BioProject [PRJNA1445501](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1445501) |
| Sample metadata | This repository, `metadata.txt` |

---

## Metadata

The sample metadata file (`metadata.txt`) is a tab-delimited text file with the following columns

| Column | Description |
|---|---|
| `sample` | Unique sample identifier |
| `fq1` | Forward FASTQ filename |
| `fq2` | Reverse FASTQ filename |
| `timepoint` | Sampling timepoint (T0 to T4 corresponding to Week 0, 2, 4, 6, 8) |
| `type` | Sample compartment (Gut, Soil, Water) |
| `shrimp` | Shrimp individual within pond |
| `pond` | Pond identifier |
| `code` | Sampling code |
| `om_percent` | Organic matter percentage |
| `nitrite` | Nitrite concentration |
| `ph` | pH |
| `ammonia` | Ammonia concentration |
| `feedinput` | Feed input |
| `temperature` | Water temperature |
| `do` | Dissolved oxygen |
| `trial` | Trial number |
| `shared` | Shared sample indicator |
| `trt` | Treatment group |

---

## Requirements

### R Environment

| Software | Version |
|---|---|
| R | 4.5.2 (2025-10-31) |
| RStudio | 2026.1.0.392 (Apple Blossom) |
| Platform | macOS Tahoe 26.3.1, aarch64-apple-darwin20 |

### R Packages

| Package | Version |
|---|---|
| `ANCOMBC` | 2.12.0 |
| `ape` | 5.8.1 |
| `Biostrings` | 2.78.0 |
| `dada2` | 1.38.0 |
| `DECIPHER` | 3.6.0 |
| `digest` | 0.6.37 |
| `emmeans` | 2.0.1 |
| `flextable` | 0.9.7 |
| `future.apply` | 1.11.3 |
| `gghalves` | 0.1.4 |
| `ggrepel` | 0.9.6 |
| `ggtext` | 0.1.2 |
| `ggVennDiagram` | 1.5.7 |
| `janitor` | 2.2.1 |
| `microbiome` | 1.32.0 |
| `multcomp` | 1.4.29 |
| `multcompView` | 0.1.10 |
| `patchwork` | 1.3.2 |
| `phangorn` | 2.12.1 |
| `phyloseq` | 1.54.2 |
| `phytools` | 2.4.4 |
| `picante` | 1.8.2 |
| `RColorBrewer` | 1.1.3 |
| `readr` | 2.1.5 |
| `scales` | 1.3.0 |
| `tidyverse` | 2.0.0 |
| `vegan` | 2.7.3 |

### External Tools

| Tool | Version | Purpose |
|---|---|---|
| VSEARCH | 2.30.0 | Chimera detection (01_16S_Pipeline.R) |
| VeryFastTree | 4.0.5 | Phylogenetic tree construction (01_16S_Pipeline.R) |

---

## How to Run

Scripts must be run in numbered order. Each script reads from and writes to an `output/` directory. Update file paths in the configuration section at the top of each script before running.

| Script | Input | Output |
|---|---|---|
| `01_16S_Pipeline.R` | `seqtab_nochim.rds` | `step4_phyloseq_object.rds`, QC summary CSVs |
| `02_Rarefaction.R` | `step4_phyloseq_object.rds` | Rarefaction plots (PNG), rarefaction CSVs |
| `03_ThresholdAnalysis.R` | `seqtab_nochim.rds` | Filtering impact matrix CSVs |
| `04_AlphaDiversity.R` | `step4_phyloseq_object.rds` | Alpha diversity plots (PNG), statistics tables (DOCX), CSVs |
| `05_BetaDiversity.R` | `step4_phyloseq_object.rds` | PCoA plots (PNG), PERMANOVA/BETADISPER CSVs |
| `06_VennDiagram.R` | `step4_phyloseq_object.rds` | Venn diagram plots (PNG) |
| `07_RelativeAbundance.R` | `step4_phyloseq_object.rds` | Alluvial plots (PNG) |
| `08_DifferentialAnalysis.R` | `step4_phyloseq_object.rds` | Volcano plots (PNG), differential abundance CSVs |

Scripts `02` and `03` can be run independently of each other after `01`. Scripts `04` through `08` all take `step4_phyloseq_object.rds` as input and can be run in any order after `01`.

---

## Contact

For questions regarding the analysis scripts, please open an issue on this repository or contact khanhnguyen@auburn.edu.
