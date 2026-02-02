This repository contains the code used to analyze data for the manuscript entitled "CSF turnover reshapes biomarker interpretation in neurodegeneration studies". Code was run on R version 4.1.1.

Scripts 01.Data_preprocessing_andPCA.R and 02.Data_analysis_ACE_cohort.R run totally or partially on omics or patient data, which can not be made publicly available due to ethical and legal reasons. 

01.Data_preprocessing_andPCA.R contains code for lipidomics and proteomics QC, processing and PCA. Furthermore, it contains code for processing metadata to generate an analysis-ready file. 
02.Data_analysis_ACE_cohort.R contains statistical analyses conducted locally using omics and metadata within the ACE CSF cohort.
03.VVIA_replications.R contains code for VVIA replication assessment using provided QTL summary statistics data
04.GNPC_data_processing_and_analysis.R contains code used within the GNPC server to process data and generate result files.
05.Result_integration.R contains code for integrating results from ACE CSF, KnightADRC and GNPC cohorts
06.AB42_and_ptau_signatures.R contains code for analyzing p-tau and AB42 CSF proteomic signatures.

The directory GWAS_scripts contains code used to analyze genomics data. To access such data, contact the corresponding author (aruiz@fundacioace.com), and Victoria Fernandez (vfernandez@fundacioace.org).

To access omics and patient data from the ACE CSF cohort (https://www.acebarcelona.org/), contact the corresponding author (aruiz@fundacioace.com) and Victoria Fernandez (vfernandez@fundacioace.org).

To access Knight ADRC CSF AD biomarker and proteomics data contact Carlos Cruchaga (Washington University, Saint Louis, MO, USA; https://knightadrc.wustl.edu).

The harmonized GNPC data has been made available for public request by the AD Data Initiative. Members of the global research community will be able to access the metadata and place a data use request via the AD Discovery Portal (https://discover.alzheimersdata.org/). Access is contingent upon adherence to the GNPC Data Use Agreement and the Publication Policies. The GNPC V1 harmonized data set (HDS) request link can be found on the GNPC website (https://www.neuroproteome.org/harmonized-data-set-hds).
