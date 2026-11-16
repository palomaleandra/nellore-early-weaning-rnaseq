# nellore-early-weaning-rnaseq

R scripts for the article:

**Skeletal muscle transcriptomic responses to early weaning in Nellore cattle: short- and long-term patterns and associations with carcass traits**
Melo PLG, Torrecilhas JA, Tinoco GLB, Russo GH, Andrini LS, Miranda LT, Figueira Filho SB, Curi RA, Moriel P, Baldassini WA, Chardulo LAL, Pereira GL.
Submitted to *Functional & Integrative Genomics*.

## Contents

| Script | Description |
|---|---|
| `R/00_setup.R` | Packages, paths, seed and analysis parameters (sourced by every script) |
| `R/01_phenotypes.R` | Performance and carcass traits (Table 1), power analysis and TOST equivalence (Tables S5 and S6) |
| `R/02_preprocess_rnaseq.R` | Low-count filtering and variance-stabilizing transformation (VST) of RNA-seq counts |
| `R/03_splsda.R` | sPLS-DA at weaning (205 d) and pre-slaughter (914 d): tuning, cross-validation and permutation tests |
| `R/04_spls_carcass.R` | Integrative sPLS between pre-slaughter transcripts and carcass traits |
| `R/06_enrichment_overlap.R` | Over-representation analysis (enrichR) and overlap between models (Euler diagram) |
| `R/07_figure_stability.R` | Selection-frequency plots of the sPLS-DA models |
| `R/run_all.R` | Runs the full pipeline in order and saves `results/sessionInfo.txt` |

## How to run

1. Open the project folder in R (the scripts use the `here` package to build paths).
2. Place the input files in `data/raw/` (`counts.tsv`) and `data/metadata/` (`sample_metadata.csv`, `phenotypes.csv`, `bodyweight_long.csv`).
3. Run `source("R/run_all.R")`. Results are written to `results/` and figures to `figures/`.

Analyses were run in R 4.3.1. Package versions are recorded in `results/sessionInfo.txt`.

## Data availability

The raw RNA-seq data (FASTQ files) are deposited in Zenodo under restricted access (https://doi.org/10.5281/zenodo.XXXXXXX). Phenotypic data are available from the corresponding author upon reasonable request. No data are included in this repository.

## Contact

Paloma Leandra Garcia Melo – paloma.leandra@unesp.br

## License

MIT
