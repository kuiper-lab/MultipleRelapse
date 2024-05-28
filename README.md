# Multiple Relapse Study

This repository contains scripts used for data analyses and visualizations in the multiple relapse study. More information about the study can be found in the [manuscript] (hyperlink follows)

## Scripts

Data analyses are conducted in Bash v4.4.20 and R v4.1.2.


### MR_signatureAnalysis_SBS.R

Script used for SBS signature analysis. Filter and cluster WGS mutations and generate a mutation matrix for de novo signature extraction. Visualise allele frequencies, mutation profiles and signature contributions.   

### MR_signatureAnalysis_DBSandInDel.R 

Script used for doublet base substitutions (DBS) and InDel signature analysis. Filter WGS mutations and generate a mutation matrix for each time point. Conduct a strict refit to calculate signature contributions.   


### denovoextraction.R

Script used for conducting non-negative matrix factorization (NMF) to perform de novo extraction. The mutation matrix of single nucleotide variants (SNVs), as created by **MR_signatureAnalysis_SBS.R**, is required as input.   


### generateDeNovoExtractionJob.sh

Wrapper script to generate jobs for de novo extraction.   


### deepSequencingOverview.R

Script used for WGS and deep sequencing data analysis. Validate WGS MAFs with deep sequencing data. Get an overview of the WGS clusters and their presentation in deep sequencing data. calculate validation rates and validate the acquired vs subclonal distinction made in the clustering of WGS mutations.   


### createdFunctions_HM-ALL.R

Script with custom functions created to parse data in the multiple relapse study.
