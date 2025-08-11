```mermaid
flowchart TB
  %% Nodes
  asm["Genome Assembly"]
  ann["Genome Annotation"]
  vc["Variant Calling"]
  vf["Variant Filtering"]
  wga["Whole Genome Alignment (Seabirds)"]
  adna["aDNA Variant Calling"]
  beast["Phylogenetics (BEAST)"]
  arg["Ancestral Recombination Graphs"]
  mig["Migration"]
  gwas["Genome-Wide Association Studies (GWAS)"]
  pgs["Population Genomic Statistics"]
  phen["Phenotypes / Health Data"]

  %% Edges
  asm --> ann
  asm --> vc --> vf
  asm --> wga
  vf --> mig
  vf --> arg
  vf --> beast
  wga --> beast
  adna --> beast
  phen --> gwas
  vf --> gwas
  vf -.-> pgs

  %% Styling
  classDef core fill:#e8f1ff,stroke:#4a90e2,stroke-width:1px,color:#0b3d91;
  classDef analysis fill:#fff4e5,stroke:#f5a623,stroke-width:1px,color:#5a3d00;
  classDef external fill:#fbe9eb,stroke:#d0021b,stroke-width:1px,color:#7a0010;

  class asm,ann,vc,vf core;
  class wga,adna,beast,arg,mig,gwas,pgs analysis;
  class phen external;
```

# Propose Directory Structure

## 'Foundational' Data Analysis (used in every analysis below)
00_Genome_Assembly/
01_Genome_Annotation/
02_Variant_Calling/
03_Variant_Filtering/

## Things that aren't used by every analysis / are their own analyses
ANALYSIS_aDNA_Variant_Calling/
ANALYSIS_Ancestral_Recombination_Graphs/
ANALYSIS_Genome_Wide_Association_Studies/
ANALYSIS_Migration/
ANALYSIS_Phylogenetics_BEAST/
ANALYSIS_Population_Genomic_Statistics/
ANALYSIS_Whole_Genome_Alignment_Seabirds/

[See Kakapo Repo](https://github.com/genomicsAotearoa/kakapo)
