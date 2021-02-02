# covid_metabolism

This repository contains the code for reproducing the analyses in:

Cheng et al. Genome-scale metabolic modeling reveals SARS-CoV-2-induced host metabolic reprogramming and identifies metabolic antiviral targets.

## 1. Prerequisites

### Softwares

* R 3.6
* R packages
  - ImNotaGit/my.utils
  - ruppinlab/Rcplex2: required to run genome-scale metabolic modeling (GEM); otherwise may be omitted
  - ruppinlab/gembox
  - other needed packages can be obtained from CRAN or Bioconductor, please see the R scripts
* IBM ILOG CPLEX Optimization Studio 12: required to run GEM; otherwise may be omitted

### Data

Download the data files from [here](url) then decompress them into the **data** folder; or it may be possible to use the dnload.sh script in the **data** folder.

## 2. Description of folders and files

The scripts should be run in the same order as they are introduced below.

### data

All required data should be saved in this folder.

* dnload.sh: a bash script for downloading the required data
* collect.validation.data.R: prepare the data for validating the MTA prediction

### expression

Scripts for gene expression-level analysis: 

* de.and.gsea.R: differential expression (DE) and gene set enrichment analysis (GSEA) between the SARS-CoV-2-infected samples and non-infected controls in each dataset 

### GEM

Scripts for genome-scale metabolic modeling (GEM) analysis:

* prepare.data.R: prepare data for GEM
* imat.and.mta.R: run iMAT and MTA on each dataset
* dflux.R: differential flux analysis between the SARS-CoV-2-infected samples and non-infected controls in each dataset 
* collect.results.R: collect the differential flux analysis and MTA results across datasets

### analyze_results

Scripts and R notebooks for inspecting the results from various analyses and for validating the MTA predictions:

* sth.R

