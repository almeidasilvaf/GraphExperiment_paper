## Overview

This repository contains code and data used in the paper
*GraphExperiment: infrastructure for high-dimensional quantitative data and
associated networks in R/Bioconductor*. Analyses are organized as a Quarto
book, with each chapter demonstrating a use case of the
[GraphExperiment](https://github.com/almeidasilvaf/GraphExperiment) package
for a different omics data type:

1. Gene coexpression and regulatory networks (bulk RNA-seq)
2. Protein co-abundance networks (single-cell proteomics)
3. Microbial co-occurrence networks (16S rRNA amplicon sequencing)
4. Cell-cell graphs for clustering (scRNA-seq)
5. Cell-cell communication networks (scRNA-seq)
6. Comparative genomics with species and phylogenetic profile graphs
7. Comparative methods on phylogenetic networks

## Repository structure

- `chapters/`: Quarto documents for each use case.
- `code/`: helper R functions sourced in the chapters.
- `data/`: small data files not available in R/Bioconductor packages.

