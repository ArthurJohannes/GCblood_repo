# GCblood


### Project description

 

A research project  to  determine signatures of glucocorticoid  (GC) regulated gene expression from public blood gene expression data.  This repository  accompanies the method section of a scientific paper entitled  **“A cortisol driven gene expression signature from circulating monocytes and neutrophils during severe inflammation”**.

Two different  codes calculate ranked averaged gene expression correlation profiles for gene queries, in multiple datasets using meta analysis. One code selects best correlation from 30 datasets, second code uses all datasets in either of 2 specific collections (each collection 15 datasets). Results from multiple gene queries with related expression (gene modules) can be combined into a single overall signature. Additionally , the project includes plotting scripts,  to show the distibution of gene (module) expression values between experimental groups,  in example datasets.  All code has been precisely annotated.

 

The meta analysis approach in the project solves problems of obtaning information for GC-regulated genes in vivo, when using single datasets, and does not require defined experimental groups for patient GC treatment.  Older transcriptomic data using microarray platforms can be efficiently reused, with automatic selection of adequate gene probes.

 

More plotting  scripts will be added to the project, for  visualizing specific gene (modules) expression in selected longitudinal experiments.  Plotting such data will be realized in a separate project using a Shiny  App. 

 

### How to install and run

 

Clone GCblood_repo (about 800MB), and run R code with provided gene expression dataset examples (30) locally in R studio environment.  R code for plotting depends on gridExtra and ggplot2. Scripts for  obtaining gene expression correlation profiles are in the code folder. Scripts for plotting are in subfolder scripts in folder plotting.

 

### How  to use project

 

Check gene expression correlation characteristics in provided blood transcriptomics collections for all your genes of interest.  Locally add  (collections of) new datasets in GEO datasets  (GDS) format.  

 

### Credits

 

Credits to the many authors of transcriptomic datasets, used here, which are mentioned in the corresponding GSE accession information at the GEO website.

 

### License



project with GPL license 



### How to contribute to project

 

Any contribution to widen such a meta analysis approach to include other types of omics datasets, or other suggestions, will be appreciated. 