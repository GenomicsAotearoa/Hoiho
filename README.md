# Hoiho

This repo contains scripts and working directories for Hoiho (yellow-eyed penguin) genomics analyses.

Authors: Joseph Guhlin<sup>1,2</sup>, Jordan Douglas<sup>3,4</sup>

# Repository layout

## Core pipeline workspaces

- `00_Genome_Assembly/` — genome assembly workspace (currently empty placeholder).
- `01_Genome_Annotation/` — genome annotation workspace (currently empty placeholder).
- `02_Variant_Calling/` — variant calling workspace (currently empty placeholder).
- `03_Variant_Filtering/` — variant filtering workspace (currently empty placeholder).

## Analysis workspaces

- `ANALYSIS_aDNA_Variant_Calling/` — ancient DNA variant calling analysis workspace (currently empty placeholder).
- `ANALYSIS_Ancestral_Recombination_Graphs/` — ARG analyses workspace (currently empty placeholder).
- `ANALYSIS_Genome_Wide_Association_Studies/` — GWAS analyses workspace (currently empty placeholder).
- `ANALYSIS_Imputation_and_Phasing/` — imputation/phasing analyses workspace (currently empty placeholder).
- `ANALYSIS_Migration/` — migration analyses workspace (currently empty placeholder).
- `ANALYSIS_Phylogenetics_BEAST/` — BEAST phylogenetics analyses workspace (currently empty placeholder).
- `ANALYSIS_Population_Genomic_Statistics/` — population genomics statistics workspace (currently empty placeholder).
- `ANALYSIS_Whole_Genome_Alignment_Seabirds/` — seabird WGA workspace (currently empty placeholder).

## Scripts and notebooks

- `scripts/` — analysis scripts, notebooks, and small utilities.
  - `scripts/general/` — exploratory notebooks and helpers (PCA, heterozygosity, sex chromosome analyses; `.ipynb`, `.sh`, `.bash`, `.qmd`, `pixi.toml`).
  - `scripts/gwas/` — phenotype processing + GWAS/TASSEL runners (`.ipynb`, `.sh`, `.py`).
  - `scripts/phylogenetics_beast/` — building sequences/SNP inputs for BEAST + utilities (`.ipynb`, `.py`, `.sh`).
    - `scripts/phylogenetics_beast/tools/` — small FASTA/region stats helpers (`.py`).
    - `scripts/phylogenetics_beast/adna/` — aDNA-specific helpers (currently just subfolders).
      - `scripts/phylogenetics_beast/adna/richdalei/` — richdalei aDNA calling/error estimation/rescue scripts (`.py`, `.sh`).
      - `scripts/phylogenetics_beast/adna/waitaha/` — waitaha aDNA calling/error estimation scripts (`.py`, `.sh`).
  - `scripts/migration/` — migration model runs + trace inspection (`.sh`, `.ipynb`, `pixi.toml`).
  - `scripts/arg/` — ARG sampling/filtering helpers (`.sh`, `.py`, `pixi.toml`).
  - `scripts/circos/` — Circos track generation notebooks + config (`.ipynb`, `circos.conf`).
  - `scripts/population_genomics/` — local PCA notebook + environment (`.ipynb`, `pixi.toml`).
  - `scripts/variant_filtering/` — standardized filtering run scripts (`.sh`).
  - `scripts/adna/` — reserved for aDNA scripts (currently empty placeholder).
  - `scripts/wga/` — whole-genome alignment helpers and utilities.
    - `scripts/wga/bigtree/` — WGA “big tree” helpers (notebook + merge script).
    - `scripts/wga/reprise/` — repeat evaluation utilities (`.py`).

# Supporting GitHub Repositories
[Modified ATLAS to reduce memory leaks and support resume functionality](https://github.com/jguhlin/ATLAS)

[See Kakapo Repo](https://github.com/genomicsAotearoa/kakapo) (will be moved to GA when final)
