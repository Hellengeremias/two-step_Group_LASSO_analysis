# Two-step Group LASSO analysis

This repository contains the R scripts used in the study:

**Selecting Genetic Variants and Interactions Associated with Amyotrophic Lateral Sclerosis: A Group LASSO Approach**

The analysis implements a two-step approach for identifying genetic variants and SNP–SNP interactions associated with amyotrophic lateral sclerosis (ALS).

## Analysis workflow

The procedure consists of two main steps:

1. **Group LASSO** is used to select candidate SNPs.
2. **glinternet** is used to evaluate main effects and pairwise interactions among the selected SNPs.

The analysis also includes bootstrap resampling and random subsets of predictors to assess the stability of selected variants and interactions.

## Repository files

* `two_step_group_lasso.R`: main script implementing the two-step analysis.
* `results_summary.Rmd`: script used to summarize and aggregate the results.

## Data availability

We analyzed data from a case-control GWAS from the National Institute of Neurological Disorders and Stroke Repository available for download through the database of Genotypes and Phenotype (dbGaP) 
Authorized Access System (dbGaP study accession phs000101.v3.p1; Data Access Request number 87433-1).

Users interested in applying the workflow to other datasets will need to adapt the data-loading and variable-mapping sections of the scripts.

## Citation

If you use or adapt this code, please cite:

Feronato SG et al. Selecting Genetic Variants and Interactions Associated with Amyotrophic Lateral Sclerosis: A Group LASSO Approach. *Journal of Personalized Medicine*. 2022;12(8):1330.

DOI: https://doi.org/10.3390/jpm12081330

Article: https://pmc.ncbi.nlm.nih.gov/articles/PMC9410070/
