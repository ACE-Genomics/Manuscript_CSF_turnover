#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args <- c("/nas/HARPONE/QTL_Analysis/pQTL/MuraliBiggs_Replication/20230531_SOMAscan_foundMurali_n19proteins.txt", "Sex_1m_2f,Age_LP")
# args <- c("/nas/HARPONE/QTL_Analysis/mQTL/02.Apply_Cruchaga_QC_extraCovars_june2023/Table_for_GWAS_Cruchaga_QC_June2023.txt", "Sex_1m_2f,Age_LP")

input.table<-args[1]
input.dir <- dirname(input.table)
covars <- unlist(strsplit(args[2], ","))

setwd(input.dir)

gwas.table <- read.table(input.table, h=T)
lipid.species <- read.table(paste(input.dir,"lipid_species.txt", sep="/"))[,1]
lipid.species <- lipid.species[lipid.species%in%names(gwas.table)]
print(paste("ANALYTES REMAINING: ", length(lipid.species)))

## Guardar covariables:

covs.table <- gwas.table[,c("IID", covars)]
print(paste("CHECK, COVARIATES ARE COMPLETE:", sum(complete.cases(covs.table)) == nrow(covs.table)))
print(paste("NUMBER OF INDIVIDUALS: ", sum(complete.cases(covs.table))))

write.table(covs.table, paste("Covs_noPCs", basename(input.table), sep="_"), row.names=F, quote=F, sep="\t")
## Guardar lipidos

for(i in 1:length(lipid.species)){
  lp.i <- lipid.species[i]
  write.table(gwas.table[!is.na(gwas.table[,lp.i]),c("IID", lp.i)], paste(input.dir, "/GWAS/", lp.i, "/", "Pheno.", lp.i, ".txt",sep=""), row.names=F, quote=F, sep="\t")
}



