# Immunity profiles in HNSCC

## Description

This code has been developed to investigate the immunity lanscape of primary head & neck squamous cell carcinomas, as presented in the manuscript "Immune cell abundance and T cell receptor landscapes suggest new patient stratification strategies in head and neck squamous cell carcinoma" (Secrier et al). 

## Code
The scripts employed in this analysis are as follows:

**immuneDeconvolution.R** - performs immune deconvolution in the cohort from expression data and defines the immunity groups

**immuneValidationInTCGA.R** - validates the uncovered immune groups in TCGA

**immunityAnalysisByGroup.R** - further analysis of the immunity groups in relation to clinical and genomic parameters

**generalImmunityAnalysis.R** - analysis of TLS enrichment in the immune groups, and links to molecular measurements via IHC

**snvs_oncoprint.Rmd** - generates an oncoprint of genomic changes in the overall cohort

**mutSignatures_analysis.Rmd** - infers mutational signatures acting in the tumours analysed in this study

**snvs_oncoprint_RTK.Rmd** - generates an oncoprint of genomic changes across genes in the Ras/ERK and PI3K/AKT pathways

**expressionAnalysis_CTME_highVSlow.Rmd** - modelling of immunity based on expression changes in genes across the Ras/ERK and PI3K/AKT pathways

**coxPH_analysis_CTME.Rmd** - Cox Proportional Hazards modelling of patient outcomes taking into account immunity

## System Requirements
Operating system(s): Unix (linux, mac)
Programming Language: R

# Copyright
This code is free and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the GNU General Public License for more details.
