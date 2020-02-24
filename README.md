
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Correlation Tree Analysis

<!-- badges: start -->

[![Codacy
Badge](https://api.codacy.com/project/badge/Grade/ba04cd22d16047bb831608b9a7a6702f)](https://www.codacy.com/app/abichat/correlationtree_analysis?utm_source=github.com&utm_medium=referral&utm_content=abichat/correlationtree_analysis&utm_campaign=Badge_Grade)
![Last-changedate](https://img.shields.io/badge/last%20change-2020--02--24-yellowgreen.svg)
[![Journal](https://img.shields.io/badge/published-bioRxiv-blue)](https://www.biorxiv.org/content/10.1101/2020.01.31.928309v1)
<!-- badges: end -->

This repository contains analysis done in the paper **Incorporating
phylogenetic information in microbiome abundance studies has no effect
on detection power and FDR control**
([bioRxiv](https://www.biorxiv.org/content/10.1101/2020.01.31.928309v1)).

Some results might be slightly different from those in the article due
to seed choice or limited number of replications in simulations.

## Structure of the repository

### Forest

Each subfolder is named after the studied dataset and contains several
files:

  - a script that performs the analysis (`analysis_dataset.R`),
  - several `.png` figures corresponding to different visualizations and
    distances,
  - saved intermediary results to avoid long recomputation (`.rds`),
  - evetually biom (`.biom`) or newick (`.nwk`) files for input data
    such as abundnace tables or phylogeny.

Each R script:

1.  loads the data,
2.  generates the forest of trees,
3.  computes pairwise distances between trees (with BHV and RF),
4.  performs PCoAs (one per distance),
5.  draws individual plots:
    1.  distance from each tree to the correlation tree,
    2.  projection of the forest on the two first axis.

### Real Datasets

This folder contains scripts to do differentialy abundance studies on
datasets Chaillou, Chlamydiae and Zeller (genus and MSP level).

As previously, it contains the R script, plots (`.png`), intermediary
results (`.rds`) and envetually biological input data (`.biom` and
`.nwk`).

Each script compare correlation and taxonomy or phylogeny in term of
detected species. It uses hierarchical FDR (from `structSSI`) for both
Chaillou and Chlamydiae datasets, and z-scores smoothing (from
`StructFDR`) for datasets from Zeller.

### Simulations

This folder contains scripts that simulates datasets according to
parametric and non parametric schemes.

The parametric simulation mimics the scheme done in Xiao et al. (2018).
It fits a negative-binomial Dirichlet-multinomial distibution on Wu
dataset and generates new differentially abundance datasets.

The non-parametric simulation uses a real dataset from Brito et al.
(2016). It generates differentially abundant species by applying
specified fold-change in half of the samples.

For each set of parameters, at least 600 replication were done. As this
is really time consumming, preprocessed data are saved in `.rds` files.

### Figures

This folder contains scripts that draw every figure in the article. Each
script takes its input in the corresponding folders of the repository.

## Reproducibility and packages

Analysis ran under R version 3.6.1.

| Package                | Version    |
| :--------------------- | :--------- |
| ape                    | 5.3        |
| biomformat             | 1.12.0     |
| broom                  | 0.5.4      |
| correlationtree        | 0.0.0.9003 |
| cowplot                | 1.0.0      |
| curatedMetagenomicData | 1.14.1     |
| distory                | 1.4.3      |
| dplyr                  | 0.8.4      |
| evabic                 | 0.0.0.9007 |
| forcats                | 0.4.0      |
| furrr                  | 0.1.0      |
| ggplot2                | 3.2.1      |
| ggstance               | 0.3.3      |
| ggtree                 | 1.16.6     |
| glue                   | 1.3.1.9000 |
| igraph                 | 1.2.4.2    |
| janitor                | 1.2.1      |
| phyloseq               | 1.28.0     |
| purrr                  | 0.3.3      |
| readr                  | 1.3.1      |
| scales                 | 1.1.0      |
| stringr                | 1.4.0      |
| StructFDR              | 1.3        |
| structSSI              | 1.1.1      |
| tidyr                  | 1.0.2      |
| tidyverse              | 1.3.0      |
| yatah                  | 0.0.1      |

To install non-CRAN packages, run these lines:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomformat")
BiocManager::install("curatedMetagenomicData")
BiocManager::install("ggtree")
BiocManager::install("phyloseq")
BiocManager::install("multtest") # dependancy for structSSI

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("abichat/correlationtree")
remotes::install_github("abichat/evabic")
remotes::install_url("https://cran.r-project.org/src/contrib/Archive/structSSI/structSSI_1.1.1.tar.gz") # archived from CRAN
```
