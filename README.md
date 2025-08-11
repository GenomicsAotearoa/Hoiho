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
