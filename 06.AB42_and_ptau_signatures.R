
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(purrr)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)
library(VennDiagram)
library(ggpubr)
library(meta)
library(ggVennDiagram)
library(forcats)

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

## Import Somascan Annots
path.annot <- "data/SOMAscan_Assay_7K_annotations_v2.xlsx"
somascan.annot <- read.xlsx(path.annot) %>% rename(Analytes=AptName) %>% select(Analytes, EntrezGeneSymbol, UniProt, Target, TargetFullName)

######## Preliminary tests  ########

## List of reproducible proteins (2395 prots)

list.good.proteins <- read.xlsx("data/SOMAscan_GoodProteins_QC.xlsx")
list.good.proteins$Analytes <- gsub("seq\\.", "seq", list.good.proteins$AptName)

## Import Muhammad Ali sAD Neuron results (2025)

mohamm.2025.results <- read.xlsx("data/AB42_ptau_signatures/Ali_2025_Neuron_AD_CSF_DAPs.xlsx",
                                 startRow = 4)[,1:31] %>% 
  rename(Analytes=Analyte) %>% mutate(Analytes=gsub("^X", "seq", Analytes))

names(mohamm.2025.results)[4:7] <- paste(c("Estimate", "SE", "P", "FDR"), "Discovery", sep="_")
names(mohamm.2025.results)[8:11] <- paste(c("Estimate", "SE", "P", "FDR"), "Replication", sep="_")
names(mohamm.2025.results)[17:20] <- paste(c("Estimate", "SE", "P", "FDR"), "Amyl.PET", sep="_")
names(mohamm.2025.results)[24:27] <- paste(c("Estimate", "SE", "P", "FDR"), "Tau.PET", sep="_")
names(mohamm.2025.results)[28:31] <- paste(c("Estimate", "SE", "P", "FDR"), "CaCo", sep="_")

## Import DIAN ADAD mutation carriers
path.cruchaga.res <- "data/AB42_ptau_signatures/mutation_carriers_DIAN.xlsx"
dian.mc.results <- read.xlsx(path.cruchaga.res)

select.mc.signif <- which(sign(dian.mc.results$Estimate.Discovery)==sign(dian.mc.results$Estimate.Replication)&dian.mc.results$FDR_p_value.Discovery<0.05&dian.mc.results$FDR_p_value.Replication<0.05)

cruchaga.mc.signif <- dian.mc.results[select.mc.signif,] %>% 
  subset(!Analytes=="*") %>%
  subset(select = -grep("Discovery|Replication", names(.))) %>% 
  mutate(Analytes=gsub("^X", "seq", Analytes))

nrow(cruchaga.mc.signif) # 240 (se pierden 2 de tau y 1 de abeta)

## APOE GNPC

gnpc.apoe <- read.xlsx("data/AB42_ptau_signatures/Supp_APOE_GNPC_paper_2025.xlsx",
                       sheet="Supplementary Table 2", startRow = 2)
names(gnpc.apoe) <- c("UniProt", "EntrezGeneSymbol", "TargetFullName", "MutualInfoValue", "Direction.APOE.GNPC")

##### IMPORT PROTEOMICS
### Import ACE proteomic models
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

ace.adbiom <- loadRData("data/AB42_ptau_signatures/AD_biomarker_models_20251001.RData")
adj.marker <- c("GAGE2A", "OPCML", "CADM2", "NPTN")

ace.adj.comparison <- lapply(adj.marker, function(adjprot){
  
  select.dfs <- c(grep("model1", names(ace.adbiom), value=T), grep(adjprot, names(ace.adbiom), value=T)) %>% sort()
  temp.obj <- ace.adbiom[select.dfs]
  
  ad.biom.models<-lapply(select.dfs, function(model){
    temp.table <- ace.adbiom[[model]] %>% rename(Estimate=Beta, Pvalue=pval, FDR=Padj_FDR) %>% 
      mutate(Analytes=gsub(".*_seq\\.", "seq", marker)) %>% 
      select(Analytes, Estimate, SE, Pvalue, FDR, zscore)
    names(temp.table)[-1] <- paste(names(temp.table)[-1], model, sep="_") %>% gsub(paste0("_", adjprot), "", .)
    temp.table
  }) %>% Reduce(function(x, y) merge(x,y,by="Analytes", all=T), .) %>% 
    merge(somascan.annot, .,by="Analytes")
  
  ## Merge APOE GNPC results
  ad.biom.models <- merge(ad.biom.models, gnpc.apoe[,c("UniProt", "MutualInfoValue", "Direction.APOE.GNPC")], by="UniProt", all.x=T)
  
  ad.biom.models
})
names(ace.adj.comparison) <- adj.marker

### Import Cruchaga proteomic models

adj.marker <- c("GAGE2A", "OPCML", "CADM2", "NPTN")
cruchaga.adj.comparison <- lapply(adj.marker, function(adjprot){
  
  input.dir <- "C:/Users/pablo/OneDrive - FUNDACIO ACE/FACE/Escritorio/PROYECTOS/HARPONE/QTLs/Analisis/06.Cruchaga_replication/2025-08_xQTL_Sumstats_FACE/2025-08_xQTL_Sumstats_FACE"
  adbiom.file.corresp <- data.frame(filename=c("CSF_Soma7K_KnightADRC_Model1_AB42_Sumstats.txt", 
                                               paste0("CSF_Soma7K_KnightADRC_Model2_AB42_Sumstats_",adjprot, ".txt"), 
                                               paste0("CSF_Soma7K_KnightADRC_Model3_AB42_Sumstats_",adjprot, ".txt"), 
                                               "CSF_Soma7K_KnightADRC_Model1_pTau_Sumstats.txt", 
                                               paste0("CSF_Soma7K_KnightADRC_Model2_pTau_Sumstats_",adjprot, ".txt"), 
                                               paste0("CSF_Soma7K_KnightADRC_Model3_pTau_Sumstats_",adjprot, ".txt")), 
                                    Tag=c("Abeta42_model1", "Abeta42_model2", "Abeta42_model3",
                                          "pTau_model1", "pTau_model2", "pTau_model3"))
  
  merged.cruchaga <- lapply(1:nrow(adbiom.file.corresp), function(i){
    fread(paste(input.dir,adbiom.file.corresp$filename[i],sep="/"),h=T, sep="\t") %>% 
      mutate(Tag=adbiom.file.corresp$Tag[i],
             zscore=Estimate/SE,
             Analytes=gsub("^X", "seq", Analytes))
  }) %>% do.call('rbind', .) %>% pivot_wider(
    id_cols = c(Analytes),
    names_from = Tag,
    values_from = c(Estimate, SE, Pvalue, FDR, zscore)
  ) %>% 
    merge(somascan.annot,., by="Analytes")
  
  ## Merge APOE GNPC results
  merged.cruchaga <- merge(merged.cruchaga, gnpc.apoe[,c("UniProt", "MutualInfoValue", "Direction.APOE.GNPC")], by="UniProt", all.x=T)
  
  merged.cruchaga[merged.cruchaga$EntrezGeneSymbol==adjprot, grep("_model2$|_model3$", names(merged.cruchaga))] <- NA
  merged.cruchaga
  
})
names(cruchaga.adj.comparison) <- adj.marker

### Check replications

### Supp table 10
lapply("OPCML", function(ref.protein){
  
  print(ref.protein)
  discovery.z.thres <- 4.496
  replication.z.thres <- 1.96
  
  discovery.data <- ace.adj.comparison
  replication.data <- cruchaga.adj.comparison
  
  temp.discovery.table <- discovery.data[[ref.protein]]
  temp.replication.table <- replication.data[[ref.protein]]
  
  lapply(c("model1", "model2", "model3"), function(model){
    
    lapply(c("Abeta42", "pTau"), function(biomarker){
      
      upreg.discovery <- temp.discovery.table$Analytes[which(temp.discovery.table[,paste("zscore", biomarker, model, sep="_")]>=discovery.z.thres)]
      downreg.discovery <- temp.discovery.table$Analytes[which(temp.discovery.table[,paste("zscore", biomarker, model, sep="_")]<=(-discovery.z.thres))]
      
      upreg.replic <- temp.replication.table$Analytes[which(temp.replication.table[,paste("zscore", biomarker, model, sep="_")]>=replication.z.thres)]
      upreg.replic <- upreg.replic[upreg.replic%in%upreg.discovery]
      
      downreg.replic <- temp.replication.table$Analytes[which(temp.replication.table[,paste("zscore", biomarker, model, sep="_")]<(-replication.z.thres))]
      downreg.replic <- downreg.replic[downreg.replic%in%downreg.discovery]
      
      temp.res <- rbind(c(length(upreg.discovery), length(upreg.replic), round(100*length(upreg.replic)/length(upreg.discovery), digits=1)),
                        c(length(downreg.discovery), length(downreg.replic), round(100*length(downreg.replic)/length(downreg.discovery), digits=1)))
      
      row.names(temp.res) <- paste(biomarker, c("Upregulated", "Downregulated"), sep="_")
      colnames(temp.res) <- paste(c("N_Disc", "N_Rep", "Perc_Rep"), model, sep="_")
      temp.res
    }) %>% do.call('rbind',.)
  }) %>% do.call('cbind',.)
  
})

## Show pheatmaps (Supp. Figure 21)

pheatmap.list <- lapply("OPCML", function(x){
  
  ace.models <- ace.adj.comparison[[x]]
  cruchaga.models <- cruchaga.adj.comparison[[x]]
  
  names(ace.models)[-2] <- paste(names(ace.models)[-1], "ACE", sep="_")
  names(cruchaga.models)[-2] <- paste(names(cruchaga.models)[-1], "KnightADRC", sep="_")
  
  cor.matrix <- merge(ace.models, cruchaga.models, by="Analytes") %>% select(matches("zscore")) %>% 
    select(matches("model1"), matches("model2"), matches("model3")) %>% 
    select(matches("Abeta"), matches("pTau")) %>% 
    cor(., use="pairwise.complete.obs")
  
  ph <- pheatmap(cor.matrix, cluster_rows = F, cluster_cols=F, display_numbers = round(cor.matrix, digits=2), 
                 fontsize_number = 10, number_color = "black", main = x,
                 color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(200), 
                 breaks = seq(-1, 1, length.out = 201),   # ensures 0 is white/neutral
  )
  ggplotify::as.ggplot(ph$gtable)
  
})

ggsave("Pheatmap_zscores_OPCML.png", width=2600, height=2200, units="px")

########  Continue with OPCML as ref protein from here   ########

ace.opcml.models <- ace.adj.comparison[["OPCML"]]
knight.opcml.models <- cruchaga.adj.comparison[["OPCML"]]

## Import and process sumstats from other studies for comparison

## Shen 2024 ADAD
path.cruchaga.res <- "data/AB42_ptau_signatures/mutation_carriers_DIAN.xlsx"
dian.mc.results <- read.xlsx(path.cruchaga.res) %>% mutate(Analytes=gsub("^X", "seq", Analytes)) %>% filter(grepl("^seq", Analytes))

select.mc.signif <- which(sign(dian.mc.results$Estimate.Discovery)==sign(dian.mc.results$Estimate.Replication)&
                            dian.mc.results$FDR_p_value.Discovery<0.05&dian.mc.results$FDR_p_value.Replication<0.05)
length(select.mc.signif) # 240, se quitan 3 (abeta42, tau y ptau)
adad.signif.analytes <- dian.mc.results$Analytes[select.mc.signif]

names(dian.mc.results)[4:7] <- paste(c("Estimate", "SE", "P", "FDR"), "Discovery", sep="_")
names(dian.mc.results)[8:11] <- paste(c("Estimate", "SE", "P", "FDR"), "Replication", sep="_")
names(dian.mc.results)[12:14] <- paste(c("Estimate", "SE", "P"), "Meta", sep="_")

analyses <- c("Discovery", "Replication", "Meta")
dian.mc.results <- lapply(analyses, function(x){
  dian.mc.results[,paste("Estimate", x, sep="_")]/dian.mc.results[,paste("SE", x, sep="_")]
}) %>% do.call('cbind', .) %>% as.data.frame() %>% setNames(paste("zscore", analyses, sep = "_")) %>% 
  cbind(dian.mc.results, .) %>% 
  mutate(Shen2025_ADAD_Direction=ifelse(Analytes%in%adad.signif.analytes&zscore_Meta>0, "Increased", 
                                        ifelse(Analytes%in%adad.signif.analytes&zscore_Meta<0, "Decreased", "No_Assoc"))) %>% 
  rename(Shen2025_ADAD_zscore_Discovery=zscore_Discovery,
         Shen2025_ADAD_zscore_Meta=zscore_Meta)

### Ali 2025, Neuron. # RPS27A- prot no detectada por Ali et al (no la tenian medida)
analyses <- c("Discovery", "Replication", "Amyl.PET", "Tau.PET", "CaCo")

mohamm.2025.results <- lapply(analyses, function(x){
  mohamm.2025.results[,paste("Estimate", x, sep="_")]/mohamm.2025.results[,paste("SE", x, sep="_")]
}) %>% do.call('cbind', .) %>% as.data.frame() %>% setNames(paste("zscore", analyses, sep = "_")) %>% 
  cbind(mohamm.2025.results, .)

## Dammer 2024. Me quedo con duplicados prot mas significativos (no reportan somamer ID). No hay dups signif en direccion contraria
path.input <- "data/AB42_ptau_signatures//Dammer2024_adn3504_Data_file_S1.xlsx"

dammer.somascan <-  read.xlsx(path.input, sheet="6.Soma DiffEx", startRow = 3)[,6:2] %>% 
  setNames(c("EntrezGeneSymbol", "logFC_CaCo", "Pvalue", "Pbonf", "Fvalue")) %>% 
  mutate(Dammer2025_Som_zscore=qnorm(Pvalue/2, lower.tail = FALSE)*sign(logFC_CaCo)) %>% 
  arrange(desc(abs(Dammer2025_Som_zscore))) %>% filter(!duplicated(EntrezGeneSymbol)) %>% # Remove duplicate Symbols, keep most significant
  filter(!EntrezGeneSymbol=="None") %>% 
  mutate(Dammer2025_Som_Direction=ifelse(Pbonf<=0.05&logFC_CaCo>0, "Increased", ifelse(Pbonf<=0.05&logFC_CaCo<0, "Decreased", "No_Assoc")))

dammer.tmt <-  read.xlsx(path.input, sheet="5.TMT-MS DiffEx", startRow = 3)[,6:1] %>% 
  setNames(c("EntrezGeneSymbol", "logFC_CaCo", "Pvalue", "Pbonf", "Fvalue", "test")) %>% 
  mutate(Dammer2025_TMT_zscore=qnorm(Pvalue/2, lower.tail = FALSE)*sign(logFC_CaCo)) %>% 
  arrange(desc(abs(Dammer2025_TMT_zscore))) %>% filter(!duplicated(EntrezGeneSymbol)) %>% # Remove duplicate Symbols, keep most significant
  filter(!EntrezGeneSymbol=="0") %>% 
  mutate(Dammer2025_TMT_Direction=ifelse(Pbonf<=0.05&logFC_CaCo>0, "Increased", ifelse(Pbonf<=0.05&logFC_CaCo<0, "Decreased", "No_Assoc")))

## Binette et al.
path.input <- "data/AB42_ptau_signatures/Binette2024_41593_2024_1737_MOESM3_ESM.xlsx"

biof2.binette <- read.xlsx(path.input, sheet="ST1_Fig2-5", startRow = 3)[,c(1,3:8)] %>% 
  setNames(c("EntrezGeneSymbol",paste(c("Estimate", "Pvalue", "FDR"), "Biofinder2_Comparison1", sep="_"),
             paste(c("Estimate", "Pvalue", "FDR"), "Biofinder2_Comparison2", sep="_"))) %>% 
  mutate(Binette_Biofinder2_zscore_Comparison1=qnorm(Pvalue_Biofinder2_Comparison1 /2, lower.tail = FALSE)*sign(Estimate_Biofinder2_Comparison1),
         Binette_Biofinder2_zscore_Comparison2=qnorm(Pvalue_Biofinder2_Comparison2 /2, lower.tail = FALSE)*sign(Estimate_Biofinder2_Comparison2)) %>% 
  mutate(Binette_Biofinder2_Comparison1_Direction=ifelse(FDR_Biofinder2_Comparison1<=0.05&Estimate_Biofinder2_Comparison1>0, "Increased", 
                                                         ifelse(FDR_Biofinder2_Comparison1<=0.05&Estimate_Biofinder2_Comparison1<0, "Decreased", "No_Assoc")),
         Binette_Biofinder2_Comparison2_Direction=ifelse(FDR_Biofinder2_Comparison2<=0.05&Estimate_Biofinder2_Comparison2>0, "Increased", 
                                                         ifelse(FDR_Biofinder2_Comparison2<=0.05&Estimate_Biofinder2_Comparison2<0, "Decreased", "No_Assoc"))) %>% 
  filter(!duplicated(.))

biof1.binette <- read.xlsx(path.input, sheet="ST2_Fig2d", startRow = 3)[,c(1,3:5)] %>% 
  setNames(c("EntrezGeneSymbol", "Estimate_Biofinder1", "Pvalue_Biofinder1", "FDR_Biofinder1")) %>% 
  mutate(Binette_Biofinder1_zscore=qnorm(Pvalue_Biofinder1/2, lower.tail = FALSE)*sign(Estimate_Biofinder1)) %>% 
  mutate(Binette_Biofinder1_Direction=ifelse(FDR_Biofinder1<=0.05&Estimate_Biofinder1>0, "Increased", 
                                             ifelse(FDR_Biofinder1<=0.05&Estimate_Biofinder1<0, "Decreased", "No_Assoc"))) %>% 
  filter(!duplicated(.))

## Tijms 2024
path.input <- "data/AB42_ptau_signatures/Tijms2024_43587_2023_550_MOESM3_ESM(1).xlsx"

tijms2024 <- read.xlsx(path.input, sheet="S Table 2", startRow = 2)[,1:15] %>% 
  select(Gene.symbol,Control, All.AD ,AD.normal.tau, AD.abnormal.tau, p.control.vs.all.AD, p.control.vs.AD.normal.tau, p.control.vs.AD.abnormal.tau) %>% 
  setNames(c("EntrezGeneSymbol", "Control", "All.AD" ,"AD.normal.tau", "AD.abnormal.tau", paste("Pvalue", c("All.AD" ,"AD.normal.tau", "AD.abnormal.tau"), sep="_")))

tijms2024[,-1] <- sapply(tijms2024[,-1], function(x){as.numeric(gsub(" .*", "", x))}) 

tijms2024 <- tijms2024 %>% select(matches("Pvalue")) %>% sapply(., p.adjust, method="fdr") %>% as.data.frame() %>% 
  setNames(paste("FDR", c("All.AD" ,"AD.normal.tau", "AD.abnormal.tau"), sep="_")) %>%
  cbind(tijms2024,.)

tijms2024 <- tijms2024 %>% 
  mutate(Tijms2024_zscore_All.AD=qnorm(Pvalue_All.AD/2, lower.tail = FALSE)*sign(All.AD),
         Tijms2024_zscore_AD.normal.tau=qnorm(Pvalue_AD.normal.tau/2, lower.tail = FALSE)*sign(AD.normal.tau),
         Tijms2024_zscore_AD.abnormal.tau=qnorm(Pvalue_AD.abnormal.tau/2, lower.tail = FALSE)*sign(AD.abnormal.tau)) %>% 
  mutate(Tijms2024_Direction_All.AD=ifelse(FDR_All.AD<=0.05&All.AD>0, "Increased", 
                                           ifelse(FDR_All.AD<=0.05&All.AD<0, "Decreased", "No assoc")),
         Tijms2024_Direction_AD.normal.tau=ifelse(FDR_AD.normal.tau<=0.05&AD.normal.tau>0, "Increased", 
                                                  ifelse(FDR_AD.normal.tau<=0.05&AD.normal.tau<0, "Decreased", "No assoc")),
         Tijms2024_Direction_AD.abnormal.tau=ifelse(FDR_AD.abnormal.tau<=0.05&AD.abnormal.tau>0, "Increased", 
                                                    ifelse(FDR_AD.abnormal.tau<=0.05&AD.abnormal.tau<0, "Decreased", "No assoc")))

## Del Campo 2022

path.input <- "data/AB42_ptau_signatures/DelCampo2022_NIHMS1902789-supplement-SourceDataTable_1.xlsx"

delcampo2022 <- read.xlsx(path.input, startRow = 5)[,1:8] %>% 
  setNames(c("UniProt", "ProtName", paste(c("Estimate", "Pvalue", "Qvalue"), "MCI.AbPlus", sep="_"), paste(c("Estimate", "Pvalue", "Qvalue"), "AD", sep="_"))) %>% 
  mutate(DelCampo_zscore_MCI.AbPlus=qnorm(Pvalue_MCI.AbPlus /2, lower.tail = FALSE)*sign(Estimate_MCI.AbPlus),
         DelCampo_zscore_AD=qnorm(Pvalue_AD /2, lower.tail = FALSE)*sign(Estimate_AD))  %>% 
  mutate(DelCampo_Direction_MCI.AbPlus=ifelse(Qvalue_MCI.AbPlus<=0.05&Estimate_MCI.AbPlus>0, "Increased", 
                                              ifelse(Qvalue_MCI.AbPlus<=0.05&Estimate_MCI.AbPlus<0, "Decreased", "No assoc")),
         DelCampo_Direction_AD=ifelse(Qvalue_AD<=0.05&Estimate_AD>0, "Increased", 
                                      ifelse(Qvalue_AD<=0.05&Estimate_AD<0, "Decreased", "No assoc")))


### Meta-analyze ACE and Knight ADRC Keep only fixed effects

# Code here:
# metaanalysis.cohorts <- lapply(c("model1", "model2", "model3"), function(model){
# 
# lapply(c("Abeta42", "pTau"), function(biomarker){
#   
#   print(biomarker)
#   
#   for.meta <- merged.opcml.results %>% 
#     column_to_rownames(var = "Analytes") %>%
#     select(matches(model)) %>% select(matches(biomarker)) %>% select(matches("Estimate"), matches("SE"))
#   
#   meta.results <- lapply(1:nrow(for.meta), function(i){
#     
#     if(i%%1000 ==0){print(i)}
#     
#     meta.summary <- metagen(for.meta[i,grep("Estimate", names(for.meta))] %>% unlist(), 
#                             for.meta[i,grep("SE", names(for.meta))] %>% unlist(), 
#                             comb.fixed = T, comb.random = T, prediction=TRUE, sm="SMD")  %>% summary()
#     
#     c(row.names(for.meta)[i],
#       meta.summary$TE.fixed,
#       meta.summary$seTE.fixed,
#       meta.summary$TE.fixed/meta.summary$seTE.fixed,
#       meta.summary$pval.fixed,
#       meta.summary$I2)
#   }) %>% do.call("rbind", .) %>% as.data.frame()
#   
#   names(meta.results) <- c("Analytes", "Estimate", "SE", "zscore", "Pvalue", "MetaHet")
#   meta.results[,-1] <- sapply(meta.results[,-1], as.numeric)
#   names(meta.results)[-1] <- paste(names(meta.results)[-1], biomarker, model, sep="_")
#   meta.results
#   
# })
# }) %>% 
#   unlist(., recursive = FALSE) %>% Reduce(function(x, y) merge(x,y,by="Analytes"), .) %>% merge(somascan.annot,., by="Analytes")
#   write.csv(metaanalysis.cohorts, "FixedEffect_Metaanalysis_ACE_KnightADRC.csv", row.names=F)
  
metaanalysis.cohorts <- read.csv("data/AB42_ptau_signatures/FixedEffect_Metaanalysis_ACE_KnightADRC.csv", h=T) %>%
    merge(., gnpc.apoe[,c("UniProt", "MutualInfoValue", "Direction.APOE.GNPC")], by="UniProt", all.x=T)

### Determine Differentially abundant proteins:

detect.daps <- function(model, input.data, z.threshold){
  
  increased.daps <- which(input.data[,paste("zscore", "pTau", model, sep="_")]>=z.threshold&input.data[,paste("zscore", "Abeta42", model, sep="_")]<=-z.threshold)
  decreased.daps <-  which(input.data[,paste("zscore", "pTau", model, sep="_")]<=-z.threshold&input.data[,paste("zscore", "Abeta42", model, sep="_")]>=z.threshold)
  
  return(list(input.data$Analytes[increased.daps], input.data$Analytes[decreased.daps]))
  
}

ace.daps.bonf <- lapply(c("model1", "model2", "model3"), detect.daps, input.data=ace.opcml.models, z.threshold=4.496) %>% purrr::flatten(.) %>% 
  setNames(c(outer(c("Increased", "Decreased"),  c("model1", "model2", "model3"), paste, sep = "_")))

ace.daps.nominal <- lapply(c("model1", "model2", "model3"), detect.daps, input.data=ace.opcml.models, z.threshold=1.96) %>% purrr::flatten(.) %>% 
  setNames(c(outer(c("Increased", "Decreased"),  c("model1", "model2", "model3"), paste, sep = "_")))

cruchaga.daps.bonf <- lapply(c("model1", "model2", "model3"), detect.daps, input.data=knight.opcml.models, z.threshold=4.496) %>% purrr::flatten(.) %>% 
  setNames(c(outer(c("Increased", "Decreased"),  c("model1", "model2", "model3"), paste, sep = "_")))

cruchaga.daps.nominal <- lapply(c("model1", "model2", "model3"), detect.daps, input.data=knight.opcml.models, z.threshold=1.96) %>% purrr::flatten(.) %>% 
  setNames(c(outer(c("Increased", "Decreased"),  c("model1", "model2", "model3"), paste, sep = "_")))

meta.analysis.daps.bonf <- lapply(c("model1", "model2", "model3"), detect.daps, input.data=metaanalysis.cohorts, z.threshold=qnorm(0.05/(2*7239), lower.tail = F)) %>% purrr::flatten(.) %>% 
  setNames(c(outer(c("Increased", "Decreased"),  c("model1", "model2", "model3"), paste, sep = "_")))

## Overlap DAPs 
# ACE: determine overlap with ADAD (DIAN cohort)
daps.dian.overlap <- lapply(c("model1", "model2", "model3"), function(model){
  
  increased.both <-intersect(ace.daps.bonf[[paste("Increased", model, sep="_")]], dian.mc.results$Analytes[dian.mc.results$Shen2025_ADAD_Direction=="Increased"])
  increased.study.only <-setdiff(ace.daps.bonf[[paste("Increased", model, sep="_")]], dian.mc.results$Analytes[dian.mc.results$Shen2025_ADAD_Direction=="Increased"])
  increased.adad.only <- setdiff(dian.mc.results$Analytes[dian.mc.results$Shen2025_ADAD_Direction=="Increased"], ace.daps.bonf[[paste("Increased", model, sep="_")]])
  
  decreased.both <-intersect(ace.daps.bonf[[paste("Decreased", model, sep="_")]], dian.mc.results$Analytes[dian.mc.results$Shen2025_ADAD_Direction=="Decreased"])
  decreased.study.only <-setdiff(ace.daps.bonf[[paste("Decreased", model, sep="_")]], dian.mc.results$Analytes[dian.mc.results$Shen2025_ADAD_Direction=="Decreased"])
  decreased.adad.only <- setdiff(dian.mc.results$Analytes[dian.mc.results$Shen2025_ADAD_Direction=="Decreased"], ace.daps.bonf[[paste("Decreased", model, sep="_")]])    
  
  list(increased.both, increased.study.only, increased.adad.only, decreased.both, decreased.study.only, decreased.adad.only)
  
}) %>% purrr::flatten(.) %>% setNames(c(outer(c("Increased_Both", "Increased_Study_only", "Increased_DIAN_only", 
                                                "Decreased_Both", "Decreased_Study_only", "Decreased_DIAN_only"), 
                                              c("model1", "model2", "model3"), paste, sep = "_")))

# Knight ADRC: overlap nominal signals with ACE significant signales
daps.knight.replication <- lapply(c("model1", "model2", "model3"), function(model){
  
  increased.replicated <-intersect(ace.daps.bonf[[paste("Increased", model, sep="_")]], cruchaga.daps.nominal[[paste("Increased", model, sep="_")]])
  increased.non.replicated <-setdiff(ace.daps.bonf[[paste("Increased", model, sep="_")]], cruchaga.daps.nominal[[paste("Increased", model, sep="_")]])
  increased.only.knight <-setdiff(cruchaga.daps.bonf[[paste("Increased", model, sep="_")]], ace.daps.bonf[[paste("Increased", model, sep="_")]])
  
  decreased.replicated <-intersect(ace.daps.bonf[[paste("Decreased", model, sep="_")]], cruchaga.daps.nominal[[paste("Decreased", model, sep="_")]])
  decreased.non.replicated <-setdiff(ace.daps.bonf[[paste("Decreased", model, sep="_")]], cruchaga.daps.nominal[[paste("Decreased", model, sep="_")]])
  decreased.only.knight <-setdiff(cruchaga.daps.bonf[[paste("Decreased", model, sep="_")]], ace.daps.bonf[[paste("Decreased", model, sep="_")]])
  
  list(increased.replicated, increased.non.replicated, increased.only.knight, decreased.replicated, decreased.non.replicated, decreased.only.knight)
  
}) %>% purrr::flatten(.) %>% setNames(c(outer(c("Increased_ACE&KnightADRC", "Increased_ACE only", "Increased_KnightADRC only", 
                                                "Decreased_ACE&KnightADRC", "Decreased_ACE only", "Decreased_KnightADRC only"), 
                                              c("model1", "model2", "model3"), paste, sep = "_")))


# Final DAPs: bonf significant in metaanalisis, nominal in both cohorts.
daps.meta.final <- lapply(c("model1", "model2", "model3"), function(model){
  
  increased.final <- intersect(ace.daps.nominal[[paste("Increased", model, sep="_")]], cruchaga.daps.nominal[[paste("Increased", model, sep="_")]]) %>% 
    intersect(., meta.analysis.daps.bonf[[paste("Increased", model, sep="_")]])
  
  decreased.final <-intersect(ace.daps.nominal[[paste("Decreased", model, sep="_")]], cruchaga.daps.nominal[[paste("Decreased", model, sep="_")]]) %>% 
    intersect(., meta.analysis.daps.bonf[[paste("Decreased", model, sep="_")]])
  
  list(increased.final, decreased.final)
  
}) %>% purrr::flatten(.) %>% setNames(c(outer(c("Increased_Meta_final", "Decreased_Meta_final"), c("model1", "model2", "model3"), paste, sep = "_")))

### Make z-score plots with no text labels 
plot.assoc.zscores.v2 <- function(input.data, models, z.threshold, list.daps){
  
  all.colors <- c("purple", "red", "blue", "darkolivegreen4", "indianred3", "navy") %>% 
    setNames(c(names(daps.dian.overlap), names(daps.knight.replication)) %>% gsub("Increased_|Decreased_", "",.) %>% gsub("_model.", "",.) %>% unique())
  
  list.plots <- lapply(models, function(model){
    
    temp.data <- input.data
    temp.daps <- list.daps[grep(model, names(list.daps))]
    
    temp.data$zscore.Abeta42 <- temp.data[,paste("zscore_Abeta42",model, sep="_")]  ## Aqui meter el cohort mas adelante
    temp.data$zscore.pTau <- temp.data[,paste("zscore_pTau",model, sep="_")]        ## Aqui meter el cohort mas adelante
    
    temp.data$Direction.APOE.GNPC <- replace_na(temp.data$Direction.APOE.GNPC, "Independent")     
    
    # Categorias por color
    color.tags <- names(temp.daps) %>% gsub("Increased_|Decreased_", "",.) %>% gsub("_model.", "",.) %>% unique()
    temp.data$Color<-NA
    for(analtye.category in color.tags){
      temp.data$Color <- ifelse(temp.data$Analytes%in%unlist(temp.daps[grep(analtye.category, names(temp.daps))]), analtye.category,temp.data$Color)
    }
    
    temp.data$Color <- factor(temp.data$Color, levels=color.tags)
    
    # color.coding <- c("purple", "red", "blue") %>% setNames(color.tags)
    color.coding <- all.colors[color.tags]
    
    data_na <- temp.data %>% filter(is.na(Color))
    data_non_na <- temp.data %>% filter(!is.na(Color))
    
    # Plot
    ggplot() +
      geom_point(data = data_na, aes(zscore.Abeta42, zscore.pTau, shape=Direction.APOE.GNPC), color = "darkgray", fill="darkgray") +
      geom_point(data = data_non_na, aes(zscore.Abeta42, zscore.pTau,color=Color, fill=Color, shape=Direction.APOE.GNPC)) +
      geom_hline(yintercept = c(-z.threshold, z.threshold), linetype = "dashed", alpha = 0.6) +
      geom_vline(xintercept = c(-z.threshold, z.threshold), linetype = "dashed", alpha = 0.6) +
      # geom_text(data = data_na, aes(x=zscore.Abeta42,y=zscore.pTau,label=EntrezGeneSymbol), color="darkgray", check_overlap =F) +
      # geom_text(data = data_non_na, aes(x=zscore.Abeta42,y=zscore.pTau,label=EntrezGeneSymbol, color=Color), check_overlap =F) +
      scale_shape_manual(
        values = c("Increased" = 24, "Decreased" = 25, "Independent" = 16),
        na.translate = TRUE
      ) +
      scale_color_manual(
        values = color.coding,
        na.translate = TRUE
      ) +
      scale_fill_manual(
        values = color.coding,
        na.translate = TRUE
      ) + theme(legend.position = "none")
  })
  
  ggarrange(plotlist=list.plots, nrow=1)
  
}

## ACE signature (Fig.6a)
plot.assoc.zscores.v2(input.data = ace.opcml.models,
                      models = c("model1", "model2", "model3"),
                      z.threshold = 4.496,
                      list.daps = daps.dian.overlap)

ggsave(filename="ADBiom_ACE_Discovery_zscorePlots_OPCML.png", width=2400, height=800, units="px")

## Knight cohort ADRC (Fig. 6b)
plot.assoc.zscores.v2(input.data = knight.opcml.models,
                      models = c("model1", "model2", "model3"),
                      z.threshold = 1.96,
                      list.daps = daps.knight.replication)
ggsave(filename="ADBiom_Cruchaga_Replication_zscorePlots_OPCML.png", width=2400, height=800, units="px")


### Make venn diagrams (using meta-analysis results)
make.venn.diagram <- function(model, input.data){
  
  protein.list <- list(input.data$Analytes[which(input.data[,paste("zscore_Abeta42", model, sep="_")] <= -4.496)],
                       input.data$Analytes[which(input.data[,paste("zscore_Abeta42", model, sep="_")] >= 4.496)],
                       input.data$Analytes[which(input.data[,paste("zscore_pTau", model, sep="_")] <= -4.496)],
                       input.data$Analytes[which(input.data[,paste("zscore_pTau", model, sep="_")] >= 4.496)])
  
  ggVennDiagram(protein.list, category.names=letters[1:length(protein.list)], label_alpha = 0.2) + 
    scale_fill_gradient(low = "white", high = "steelblue") + 
    scale_color_manual(values = rep("black", 4)) + 
    scale_alpha_manual(values=rep(0, length(protein.list))) + 
    theme(legend.position = "none") +
    ggtitle(model)
  ggsave(paste0("ggVenn_ADcorebiom_",model,".png"), height = 1000, width=1000, units="px")
  
}

lapply(c("model1", "model2", "model3"), make.venn.diagram, input.data=metaanalysis.cohorts)

###### Integrate our results and external studies. Save table for supp. #####

final.daps <- lapply(grep("model2", names(daps.meta.final), value=T), function(x){
  df <- data.frame(Analytes=daps.meta.final[[x]],Study_result=gsub("_.*", "", x))
}) %>% do.call('rbind', .) 

joined.results <- ace.opcml.models %>% select(Analytes, matches("model"), -matches("SE"), -matches("FDR")) %>% rename_with(~ paste0(.x, "_ACE"), -Analytes) %>% 
  full_join(.,knight.opcml.models %>% select(Analytes, matches("model"), -matches("SE"), -matches("FDR")) %>% rename_with(~ paste0(.x, "_Knight"), -Analytes)) %>%
  full_join(., metaanalysis.cohorts %>% select(Analytes, matches("model"), -matches("SE"), -matches("FDR")) %>% rename_with(~ paste0(.x, "_Meta"), -Analytes)) %>% 
  select(Analytes, matches("Estimate"), matches("Pvalue"), matches("zscore")) %>%
  select(Analytes, matches("_ACE"), matches("_Knight"), matches("_Meta")) %>%
  select(Analytes, matches("Abeta42"), matches("pTau")) %>%
  select(Analytes, matches("model1"), matches("model2"), matches("model3")) %>% right_join(final.daps, ., by="Analytes") %>% 
  right_join(somascan.annot,.,by="Analytes")


## Annotate APOE effect
joined.results <- merge(joined.results, gnpc.apoe[,c("UniProt", "Direction.APOE.GNPC")], by="UniProt", all.x=T)

### Annotate effect dirs
## ADAD shen 2025
joined.results <- dian.mc.results %>% 
  select(Analytes, Shen2025_ADAD_Direction, Shen2025_ADAD_zscore_Discovery, Shen2025_ADAD_zscore_Meta) %>%
  merge(joined.results,.,by="Analytes", all.x=T)

## sAD Ali 2025
joined.results <- mohamm.2025.results %>% 
  mutate(Ali2025_sAD_Direction=ifelse(zscore_Discovery>0&!is.na(weightedZ), "Increased", ifelse(zscore_Discovery<0&!is.na(weightedZ), "Decreased", "No_Assoc"))) %>% 
  rename(Ali2025_sAD_zscore_Discovery=zscore_Discovery) %>% 
  select(Analytes, Ali2025_sAD_Direction, Ali2025_sAD_zscore_Discovery) %>% 
  merge(joined.results,.,by="Analytes", all.x=T)

## sAD Dammer 2024
joined.results <- merge(joined.results, 
                        dammer.somascan[,c("EntrezGeneSymbol", "Dammer2025_Som_zscore", "Dammer2025_Som_Direction")],all.x=T, by="EntrezGeneSymbol") %>% 
  merge(., dammer.tmt[,c("EntrezGeneSymbol", "Dammer2025_TMT_zscore", "Dammer2025_TMT_Direction")], all.x=T, by="EntrezGeneSymbol")


## Binette 2024, BioFinder1 and Biofinder2
joined.results <-  merge(joined.results,biof2.binette %>% select(EntrezGeneSymbol, matches("zscore"), matches("Direction")),by="EntrezGeneSymbol", all.x=T)
joined.results <-  merge(joined.results,biof1.binette %>% select(EntrezGeneSymbol, matches("zscore"), matches("Direction")),by="EntrezGeneSymbol", all.x=T)

## Tijms 2024
joined.results <-  merge(joined.results,tijms2024 %>% select(EntrezGeneSymbol, matches("zscore"), matches("Direction")),by="EntrezGeneSymbol", all.x=T)

## Del Campo 2022
joined.results <-  merge(joined.results,delcampo2022 %>% select(UniProt, matches("zscore"), matches("Direction")),by="UniProt", all.x=T)

## Export final report. (Supp. Table 11).
write.xlsx(joined.results %>% arrange(zscore_Abeta42_model2_Meta), "Merged_Ab_Tau_signatures_forSupp.xlsx")

##### Make between-study Venn Diagrams, export for ORA  #######

### Get protein lists
quick.entrez <- function(x){return(somascan.annot$EntrezGeneSymbol[somascan.annot$Analytes%in%x] %>% unique() %>% unlist())}
quick.zscores <- function(x){print(ace.opcml.models[ace.opcml.models$Analytes%in%x,] %>% 
                                     select(EntrezGeneSymbol, matches("zscore"), Direction.APOE.GNPC, ADAD_Direction, sAD_Direction) %>% 
                                     arrange(zscore_Abeta42_model2 ))}


## AD increased DAPs (Fig 6d)

ggVennDiagram(list(joined.results$Analytes[which(joined.results$Ali2025_sAD_Direction=="Increased")] %>% quick.entrez(), 
                   joined.results$Analytes[which(joined.results$Shen2025_ADAD_Direction=="Increased")] %>% quick.entrez(), 
                   joined.results$Analytes[which(joined.results$Study_result=="Increased")] %>% quick.entrez()),
              category.names=c("sAD", "ADAD", "m2"), label_alpha = 0.2) + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  scale_color_manual(values = rep("black", 4)) + 
  scale_alpha_manual(values=rep(0, 3)) + 
  theme(legend.position = "none")
# ggsave("Venn_diagram_Increased_DAPs_Shen_Ali.png", height=800, width=800, units="px")


## AD decreased DAPs (Fig 6e)
ggVennDiagram(list(joined.results$Analytes[which(joined.results$Ali2025_sAD_Direction=="Decreased")] %>% quick.entrez(), 
                   joined.results$Analytes[which(joined.results$Shen2025_ADAD_Direction=="Decreased")] %>% quick.entrez(), 
                   joined.results$Analytes[which(joined.results$Study_result=="Decreased")] %>% quick.entrez()),
              category.names=c("sAD", "ADAD", "ACE"), label_alpha = 0.2) + 
  scale_fill_gradient(low = "white", high = "darkorange3") + 
  scale_color_manual(values = rep("black", 4)) + 
  scale_alpha_manual(values=rep(0, 3)) + 
  theme(legend.position = "none") 
# ggsave("Venn_diagram_Decreased_DAPs_Shen_Ali.png", height=800, width=800, units="px")


## Visualize enrichments (Fig. 6c)

# Export for WebGestalt

# write.table(quick.entrez(daps.meta.final[["Increased_Meta_final_model2"]]), "DAPs_Increased_finalMeta.txt", row.names=F, col.names = F, quote=F)
# write.table(quick.entrez(daps.meta.final[["Decreased_Meta_final_model2"]]), "DAPs_Decreased_finalMeta.txt", row.names=F, col.names = F, quote=F)

path.input <- "data/AB42_ptau_signatures/enrichment_results_wg_result1760096842.txt"
top.10.increased.webgestalt <- read.table(path.input, h=T, sep="\t")

num_cols <- c("size","overlap","expect","enrichmentRatio","pValue","FDR")
df_plot <- top.10.increased.webgestalt %>%
  mutate(across(all_of(num_cols), ~ as.numeric(gsub(",", ".", .x)))) %>%
  mutate(
    direction = if_else(enrichmentRatio >= 1, "positive", "negative"),
    signif = -log10(FDR)  # larger = more significant
  )

ggplot(
  df_plot,
  aes(x = enrichmentRatio, y = fct_reorder(description, enrichmentRatio))
) +
  geom_point(
    aes(size = signif, fill = direction, alpha = signif),
    shape = 21, color = "black"
  ) +
  scale_fill_manual(
    values = c(negative = "#E67E22", positive = "#2C7FB8"),
    name = "Direction"
  ) +
  # FDR-derived aesthetics: include 0 on the scale
  scale_size_continuous(limits = c(0, NA), range = c(2.5, 10),
                        name = expression(-log[10](FDR))) +
  scale_alpha_continuous(limits = c(0, NA), range = c(0.3, 1), guide = "none") +
  # X axis (enrichment ratio): include 0 on the axis
  scale_x_continuous(limits = c(1,7), expand = expansion(mult = c(0, 0.05))) +
  theme_classic(base_size = 13) +
  labs(x = "Enrichment ratio", y = "Term (GO/KEGG)")

# ggsave("DotPlot_Webgestalt_Increased_finalMeta.png", width=2200, height=1400, units="px")


###### Check Top PC1 protein z-scores (Fig. 6f)#######
# Import ACE PC1 assoc
ace.pc1.signature <- read.xlsx("data/Proteogenomic_signature.xlsx") %>% 
  mutate(Analytes=gsub("^seq\\.", "seq", AptName)) %>% 
  select(Analytes, EntrezGeneSymbol, Beta_Random.META_PC1, Beta_Random.META_PC2)

# Define studies and variables
test.studies <- data.frame(VarName=c("Shen2025_ADAD_zscore_Discovery", "Ali2025_sAD_zscore_Discovery", "Dammer2025_Som_zscore", 
                                     "Binette_Biofinder2_zscore_Comparison1", "Binette_Biofinder2_zscore_Comparison2", "Binette_Biofinder1_zscore",
                                     "Tijms2024_zscore_AD.normal.tau", "Tijms2024_zscore_AD.abnormal.tau", "DelCampo_zscore_MCI.AbPlus", "DelCampo_zscore_AD"),
                           Study=c("ADADmc Shen", "sAD Ali", "sAD Dammer", "sAD Binette A+T-", "BioFindr2_A+T-vA+T+", "BioFindr1_A-vA+", "sAD Tijms A+T-", "sAD Tijms A+T+",
                                   "sAD_DelCampo_MCIAb-", "sAD DelCampo A+T+"),
                           ID_Type=c("Analytes", "Analytes", "EntrezGeneSymbol", "EntrezGeneSymbol", "EntrezGeneSymbol", "EntrezGeneSymbol", "EntrezGeneSymbol", "EntrezGeneSymbol",
                                     "EntrezGeneSymbol", "EntrezGeneSymbol"),
                           Keep=c("yes", "yes", "yes", "yes", "no", "no", "yes", "yes", "no", "yes"))


input.data <- joined.results
extracted.zscores <- lapply(1:nrow(test.studies), function(i){
  
  if(test.studies$Keep[i]=="no"){return(NULL)}
  print(i)
  
  z.score.var <- test.studies$VarName[i]
  study.name <- test.studies$Study[i]
  id.type <- test.studies$ID_Type[i]
  
  temp.table <- merge(ace.pc1.signature, input.data[,c("Analytes", "EntrezGeneSymbol", z.score.var)], by=id.type) %>% 
    filter(complete.cases(.)) %>% arrange(desc(Beta_Random.META_PC1)) %>% mutate(Study=study.name) 
  
  names(temp.table)[names(temp.table)==z.score.var] <- "zscore"
  names(temp.table) <-gsub("\\.x$", "", names(temp.table))
  temp.table <- temp.table[!duplicated(temp.table[,id.type]),]
  
  temp.table[,c("Analytes", "EntrezGeneSymbol", "Beta_Random.META_PC1", "zscore", "Study")] 
  
}) %>% do.call('rbind', .) %>% mutate(Study=factor(Study, levels=test.studies$Study))

ggplot(extracted.zscores, aes(Beta_Random.META_PC1, zscore, color=Study, fill = Study)) + geom_smooth(alpha = 0.25, method = "gam") + 
  theme_classic() + geom_hline(yintercept=0, linetype = "dashed", alpha = 0.6) + 
  labs(x = "Beta PC1",
       y = "AD z-score")

# ggsave(filename="Trends_ADassoc_byPC1_beta.png", width=1800, height=750, units="px")

