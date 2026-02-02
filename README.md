# **CSF Turnover and Biomarker Interpretation in Neurodegeneration**

This repository contains the code used to analyze data for the manuscript entitled: “CSF turnover reshapes biomarker interpretation in neurodegeneration studies”

All analyses were performed using R version 4.1.1.

## **Repository Overview**

This repository is organized into analysis scripts corresponding to different cohorts, data types, and analytical steps. Due to ethical and legal constraints, script 01 and 02 operate on omics or patient-level data that cannot be made publicly available.

Analysis Scripts

### 01. Data preprocessing and PCA

Script: **01.Data_preprocessing_andPCA.R**

-Lipidomics and proteomics quality control (QC)

-Data processing and normalization

-Principal Component Analysis (PCA)

-Metadata processing to generate an analysis-ready dataset

Note: This script runs totally or partially on omics and patient data that cannot be shared publicly.

### 02. ACE CSF cohort data analysis

Script: **02.Data_analysis_ACE_cohort.R**

Statistical analyses conducted locally

Uses omics and clinical metadata from the ACE CSF cohort

Note: This script contains analyses based on restricted data.

### 03. VVIA replication
Script: **03.VVIA_replications.R**

Replication analyses using provided QTL summary statistics

Does not require access to individual-level patient data

### 04. GNPC data processing
Script: **04.GNPC_data_processing_and_analysis.R**

Code executed within the GNPC server environment

Processes GNPC data and generates result files

### 05. Cross-cohort result integration
Script: **05.Result_integration.R**

Integration of results from: ACE CSF cohort, Knight ADRC cohort and GNPC cohort

### 06. Aβ42 and p-tau proteomic signatures
File: **06.AB42_and_ptau_signatures.R**

Analysis of CSF proteomic signatures associated with Aβ42 and p-tau


## GWAS Analyses

The directory "GWAS_scripts" contains code used to analyze genomic data.

Access to genomic data is restricted. To request access, please contact:

Agustin Ruiz (corresponding author): aruiz@fundacioace.com

Victoria Fernandez: vfernandez@fundacioace.org

## Data Access Information

### ACE CSF cohort
Institution: Fundació ACE (https://www.acebarcelona.org/
)

Access requests:

aruiz@fundacioace.com

vfernandez@fundacioace.org

### Knight ADRC cohort
Institution: Washington University in St. Louis
Website: https://knightadrc.wustl.edu

Contact:

Carlos Cruchaga

### GNPC harmonized data
The GNPC harmonized data have been made available for public request by the AD Data Initiative.

Members of the global research community may access metadata and submit a data use request via the AD Discovery Portal:
https://discover.alzheimersdata.org/

Access is contingent upon adherence to the GNPC Data Use Agreement and Publication Policies.

The GNPC V1 Harmonized Data Set (HDS) request link is available at:
https://www.neuroproteome.org/harmonized-data-set-hds

## Notes

This repository is intended to support transparency and reproducibility of analytical workflows. Scripts that require controlled-access data are provided for methodological reference only and cannot be executed without appropriate data access.
