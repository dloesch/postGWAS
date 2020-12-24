# postGWAS
 scripts for plotting and evaluating GWAS summary statistics

 ## Plot
 R scripts with argument handling for generating GWAS plots (QQ plot, Manhattan plot, locuszoom plot)
 
 ### plot_locuszoom.R: R script that uses GWAS summary statistics and a provided genotype file to estimate LD and create a locsuzoom-style plot. 
 1. Requirements: BCFTools and PLINK installed and in path; R packages ggplot2, data.table, ggrepel, cowplot, argparser, tools; and a file defining gene boundaries (provided)
 2. Basic usage: Rscript plot_locuszoom.R -s _sumstats file_ -g _genotype file_ (See --help for full options) 
 3. Important note: This script was written to create a locuszoom-style plot using the same genotype file with which the GWAS summmary statistics were generated.

 ### plot_GWAS.R: R script for generating a Manhattan plot and QQ plot
 1. Requirements: R packages data.table, argparser, tools
 2. Basic usage: Rscript plot_GWAS.R --results _sumstats file_ (See --help for full options)

 ## Finemap 
 Simple pipeline for running PAINTOR 3.0 and the annotations compiled by the developers
 1. Requirements: BCFtools and PAINTOR 3.0 (note: run in python2 environment for best results)
 2. Workflow: (A) processes summary statistics using provided PAINTOR 3.0 utilities, (B) estiamtes base model without annotations, (C) extract only desired annotations, (D) test each annotation (stopping if one reaches statistical significance with Bonferroni correction), keeps top 3 annotations (checking for correlation), and extracts the 95% credible set. 
 3. Imporant note: scripts are formatted for my system. Modifications will be needed prior to implementing locally. 
 4. PAINTOR citation: Kichaev G, Yang W-Y, Lindstrom S, et al. Integrating Functional Data to Prioritize Causal Variants in Statistical Fine-Mapping Studies. PLOS Genetics. 2014;10(10):e1004722.



