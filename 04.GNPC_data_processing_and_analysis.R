library(openxlsx)
library(ggplot2)
library(DT)
library(dplyr)
library(kableExtra)
library(data.table)
library(dplyr)
library(pheatmap)
library(tidyr)
library(meta)
library(grid)

hed <- function(x) x[1:6, 1:6]

### 1. Import files

## my files
annot <- read.xlsx("~/files/ACE/Ref_files/20230502_SomaScan_Annotation_data_Human_Protein_Apt.xlsx", detectDates = T)
good.proteins <-  read.xlsx("~/files/ACE/Ref_files/20240429_SOMAscan_Proteomic_GoodProteins_n2395_FailQC_n33.xlsx", sheet="SOMAscan Good proteins QC WAshU", startRow = 2)

## Import proteomics data

somascan.gnpc.v1 <- fread("~/files/GNPC_HDS_general/V1.3MS/SomalogicProteomics_all.csv", h=T)
metadata.gnpc.v1 <- read.csv("~/files/GNPC_HDS_general/V1.3MS/SomalogicMetaV1_3ms.csv", h=T)
clinicaldata.gnpc.v1 <- read.csv("~/files/GNPC_HDS_general/V1.3MS/ClinicalV1_3ms.csv", h=T)
person.mapping <- read.csv("~/files/GNPC_HDS_general/V1.3MS/PersonMappingV1_3ms.csv", h=T)

metadata.gnpc.v1$Dataset <- paste0("Cont", metadata.gnpc.v1$contributor_code, "_", gsub(" ", "", metadata.gnpc.v1$sample_matrix))

## Preprocess
somascan.gnpc.v1 <- as.data.frame(somascan.gnpc.v1)
row.names(somascan.gnpc.v1) <- somascan.gnpc.v1$sample_id

test.datasets <- unique(metadata.gnpc.v1$Dataset) %>% grep("-1",. , invert = T, value=T) %>% sort() ## Define datasets to test
# test.datasets <- grep("CSF", test.datasets, value=T) ## Solo CSF

### 2. Rm non-informative, log transform, scale and remove outliers

preprocessed.gnpc <- lapply(test.datasets, function(x){
  
  temp.table <- somascan.gnpc.v1[metadata.gnpc.v1$sample_id[which(metadata.gnpc.v1$Dataset == x)],]
  protein.vals <- temp.table[,-c(1:4)]
  
  print(paste0(x, "N = ", nrow(protein.vals)))
  
  ## Remove non-informative:
  
  print(paste0(x, ", removed non informative (IQR=0) somamers: ", sum(sapply(protein.vals, IQR, na.rm=T)==0)))
  protein.vals <- protein.vals[,!sapply(protein.vals, IQR, na.rm=T)==0]
  
  ## Log10 and scale
  min.vals <- sapply(protein.vals, min)
  protein.vals.transf <- lapply(1:ncol(protein.vals), function(i){log10(protein.vals[,i]-min.vals[i]+abs(min.vals[i]))}) %>% do.call("cbind",.) %>% scale()
  
  ## Outlier removal:
  medians <- apply(protein.vals.transf, 2, median, na.rm = TRUE)
  iqrs <- apply(protein.vals.transf, 2, IQR, na.rm = TRUE)
  
  low.outlier.thres <- medians - 1.5 * iqrs
  high.outlier.thres <- medians + 1.5 * iqrs 
  
  outliers <- sweep(protein.vals.transf, 2, low.outlier.thres, `<`) |
    sweep(protein.vals.transf, 2, high.outlier.thres, `>`)
  
  print(paste0(x, ", removed outliers: ", sum(outliers), ", % values = ", round(100*sum(outliers)/(nrow(protein.vals.transf)*ncol(protein.vals.transf)), digits=1)))
  
  protein.vals.filtered <- protein.vals.transf
  protein.vals.filtered[outliers] <- NA
  colnames(protein.vals.filtered) <- names(protein.vals)
  cbind(temp.table[,c(1:4)], protein.vals.filtered)
  
})
names(preprocessed.gnpc) <- test.datasets
# save(preprocessed.gnpc, file="GNPC_preprocessed_20250804.RData")

### 3. RUN PCA (in CSF keep only reproducible proteins)

pcs.datasets <- lapply(1:length(preprocessed.gnpc), function(x){
  
  temp.table <- preprocessed.gnpc[[x]][,-c(1:4)]
  
  ## CSF keep only good prots.
  if(grepl("CSF", names(preprocessed.gnpc)[x])){temp.table <- temp.table[,names(temp.table)%in%gsub("\\.", "_", good.proteins$AptName)]}
  
  ## Calculate mean standardized values
  msv <- temp.table %>% t() %>% apply(.,2, mean, na.rm=T)
  
  print(paste0(names(preprocessed.gnpc)[x], "; PCs calculated for ", nrow(temp.table), " samples, ", ncol(temp.table), " proteins"))
  
  # Transpose and re-scale
  temp.table <- temp.table %>% t() %>% scale() 
  corr_matrix <- cov(temp.table, use="pairwise.complete.obs")
  # PCA
  data.pca <- princomp(corr_matrix)
  
  merged.pcs <- cbind(msv, data.pca$scores[,1:50])
  
  # Flip PCs (higher value, higher MSV)
  flip.pcs <- lapply(paste0("Comp.", 1:50), function(pc){coef(lm(merged.pcs[,"msv"] ~ merged.pcs[,pc]))[2]}) %>% unlist() %>% sign()
  merged.pcs[,-1] <- sweep(merged.pcs[,-1], 2, flip.pcs, `*`)  
  merged.pcs
})

names(pcs.datasets) <- test.datasets
# save(pcs.datasets, file="PCs_by_dataset_20250804.RData")

## percentage Variance explained:
lapply(1:length(pcs.datasets), function(x){  
  data.pca <- pcs.datasets[[x]]
  data.pca <- data.pca[,grep("Comp", colnames(data.pca))]
  pr.variance <- apply(data.pca, 2, sd)^2
  pve <- pr.variance/sum(pr.variance)
  plot(1:15, 100*pve[1:15], type="b", xlab="Principal component", ylab="Variance explained (%)", main=names(pcs.datasets)[x])
  
})



### 4. RUN raw linear models

# Define function
multi.lm.raw <- function(test.pheno, out.file){
  
  raw.lm.results <- lapply(test.datasets, function(x){
    
    pc.data <- pcs.datasets[[x]] %>% scale()
    temp.table <- preprocessed.gnpc[[x]][,-c(1:4)] %>% scale() %>% as.data.frame()
    
    print(paste0(x," Running somamer vs ",test.pheno," lm; ", nrow(temp.table), " samples, ", ncol(temp.table), " proteins"))
    
    results.lm <-  lapply(temp.table, function(prot){summary(lm(pc.data[,test.pheno] ~ prot))$coefficients[2,]}) %>% 
      do.call("rbind",.) %>% apply(.,2,as.numeric) %>%
      as.data.frame() %>% mutate(marker = colnames(temp.table))
    names(results.lm) <- c("Beta", "SE", "zscore", "pval", "marker")
    
    results.lm[,c("marker", "Beta", "SE", "zscore", "pval")]
  })
  names(raw.lm.results) <- test.datasets
  save(raw.lm.results, file=out.file)
  
}

## Run models

load("PCs_by_dataset_20250804.RData")
load("GNPC_preprocessed_20250804.RData")

multi.lm.raw(test.pheno="msv", out.file="Raw_lm_msv_results_20250804.RData") # MSV
multi.lm.raw(test.pheno="Comp.1", out.file="Raw_lm_Comp.1_results_20250804.RData") # PC1
multi.lm.raw(test.pheno="Comp.2", out.file="Raw_lm_Comp.2_results_20250804.RData") # PC2



### 5. METAANALYZE Plasma and CSF

## Import tables 

metadata.gnpc.v1 <- read.csv("~/files/GNPC_HDS_general/V1.3MS/SomalogicMetaV1_3ms.csv", h=T)
metadata.gnpc.v1$Dataset <- paste0("Cont", metadata.gnpc.v1$contributor_code, "_", gsub(" ", "", metadata.gnpc.v1$sample_matrix))

test.datasets <- unique(metadata.gnpc.v1$Dataset) %>% grep("-1",. , invert = T, value=T) %>% sort() ## Define datasets to test

annot.somascan <- read.csv("~/files/GNPC_HDS_general/V1.3MS/SomalogicAnalyteInfoV1_3ms.csv", h=T) # Load Annot
names(annot.somascan)[names(annot.somascan)=="column_name"] <- "marker"

## Define function
run.metaanalysis <- function(input.obj, data.selection, out.file="output.txt"){
  
  table.for.meta <- lapply(data.selection, function(x){input.obj[[x]] %>% mutate(Cohort=x)}) %>% do.call('rbind',.)
  
  meta.results<-lapply(unique(table.for.meta$marker), function(prot){
    temp.meta <- table.for.meta[table.for.meta$marker==prot,]
    m1 <- metagen(Beta, SE, data=temp.meta, common = T, random = T, prediction=TRUE, sm="SMD")
    meta.summary <- summary(m1)
    c(prot,
      meta.summary$TE.fixed,
      meta.summary$seTE.fixed,
      meta.summary$pval.fixed,
      meta.summary$TE.random,
      meta.summary$seTE.random,
      meta.summary$pval.random,
      meta.summary$I2)
  }) %>% do.call("rbind", .) %>% as.data.frame()
  names(meta.results) <- c("marker", "Beta_Fixed", "SE_Fixed", "pval_Fixed","Beta_Random", "SE_Random", "pval_Random", "I2")
  meta.results[,-1] <- sapply(meta.results[,-1], as.numeric)
  
  meta.report <- lapply(data.selection, function(x){
    temp.table <- input.obj[[x]]
    names(temp.table)[-1] <- paste(names(temp.table)[-1],x,sep=".")
    temp.table
  }) %>% Reduce(function(x, y) merge(x,y,by="marker", all=T), .) %>% 
    left_join(.,meta.results, by="marker") %>%
    left_join(annot.somascan,.,by="marker") %>% 
    mutate(zscore_Fixed=Beta_Fixed/SE_Fixed,
           zscore_Random=Beta_Random/SE_Random)
  
  write.table(meta.report, out.file)
  
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## CSF
## MSV
raw.msv.lm.results <- loadRData("Raw_lm_msv_results_20250804.RData")  

run.metaanalysis(input.obj=raw.msv.lm.results,
                 data.selection=grep("CSF", test.datasets, value=T),
                 out.file="Meta_Report_CSF_MSV.txt")

## PC1:
raw.comp.1.lm.results <- loadRData("Raw_lm_Comp.1_results_20250804.RData")  

run.metaanalysis(input.obj=raw.comp.1.lm.results,
                 data.selection=grep("CSF", test.datasets, value=T),
                 out.file="Meta_Report_CSF_PC1.txt")


## PC2
raw.comp.2.lm.results <- loadRData("Raw_lm_Comp.2_results_20250804.RData") 

run.metaanalysis(input.obj=raw.comp.2.lm.results,
                 data.selection=grep("CSF", test.datasets, value=T),
                 out.file="Meta_Report_CSF_PC2.txt")


## EDTAPlasma
## MSV
load("Raw_lm_MSV_results_20250804.RData")

run.metaanalysis(input.obj=raw.msv.lm.results,
                 data.selection=grep("EDTAPlasma", test.datasets, value=T),
                 out.file="Meta_Report_EDTAPlasma_MSV.txt")

## PC1:
load("Raw_lm_Comp.1_results_20250804.RData")
run.metaanalysis(input.obj=raw.comp.1.lm.results,
                 data.selection=grep("EDTAPlasma", test.datasets, value=T),
                 out.file="Meta_Report_EDTAPlasma_PC1.txt")


## PC2
load("Raw_lm_Comp.2_results_20250804.RData")
run.metaanalysis(input.obj=raw.comp.2.lm.results,
                 data.selection=grep("EDTAPlasma", test.datasets, value=T),
                 out.file="Meta_Report_EDTAPlasma_PC2.txt")



### Report top hits, make forest plots, etc...

make.forest.plots <- function(input.obj, data.selection, select.prot){

  table.for.meta <- lapply(data.selection, function(x){input.obj[[x]] %>% mutate(Cohort=x)}) %>% do.call('rbind',.)
  
  temp.meta <- table.for.meta[table.for.meta$marker==select.prot,]
  m1 <- metagen(Beta, SE, data=temp.meta, studlab=Cohort, sm="Beta")
  
  forest(m1,col.common="darkgoldenrod2",col.random="brown",col.diamond.common="darkgoldenrod2",col.diamond.random="brown",
         col.diamond.lines.random="brown",col.diamond.lines.common="darkgoldenrod2", overall = T, overall.hetstat = T, test.overall = T, 
         leftcols = c("studlab"), leftlabs = c("Cohort"), addrow = T, addrow.overall = T, colgap.forest.left = "15mm", 
         col.square = "cornflowerblue", print.tau2 = F, print.Q = T, backtransf = TRUE,scientific.pval = T, addrows.below.overall = 3,
         main="test")
  # grid.text(paste(annot.somascan$entrez_gene_symbol[annot.somascan$marker==select.prot], select.prot), .5, .9, gp=gpar(cex=1.5))
  print(paste(annot.somascan$entrez_gene_symbol[annot.somascan$marker==select.prot], select.prot))
}


## ggplot regressions
make.regression.plots <- function(select.pc, data.selection, select.prot){
  
  for.plot <- lapply(data.selection, function(x){
    
    if (!select.prot %in% names(preprocessed.gnpc[[x]])) return(NULL)
    temp.df <- data.frame(PC=pcs.datasets[[x]][,select.pc] %>% scale(),
                          Protein=preprocessed.gnpc[[x]][,select.prot] %>% scale(),
                          Cohort=x)
  }) %>% do.call('rbind',.)
  
  
  ggplot(for.plot, aes(PC, Protein, col=Cohort)) +
    geom_point(alpha = 0.1) +  # optional: to show the raw data
    geom_smooth(method = "lm", se = TRUE) +  # se = TRUE plots 95% CI
    # theme_minimal() +
    labs(title = paste(annot.somascan$entrez_gene_symbol[annot.somascan$marker==select.prot], select.prot),
         x = select.pc, y = select.prot)
  
}

csf.msv.meta <- read.table("Meta_Report_CSF_MSV.txt", h=T)
csf.pc1.meta <- read.table("Meta_Report_CSF_PC1.txt", h=T)
csf.pc2.meta <- read.table("Meta_Report_CSF_PC2.txt", h=T)

plasma.msv.meta <- read.table("Meta_Report_EDTAPlasma_MSV.txt", h=T)
plasma.pc1.meta <- read.table("Meta_Report_EDTAPlasma_PC1.txt", h=T)


## Test top hits

head(csf.pc1.meta$marker[order(csf.pc1.meta$Beta_Fixed, decreasing=T)])
head(csf.pc1.meta$marker[order(csf.pc2.meta$Beta_Fixed, decreasing=T)])
head(csf.pc1.meta$marker[order(plasma.pc1.meta$Beta_Fixed, decreasing=T)])



cor.matrix <- csf.pc1.meta %>% select(matches("Beta.Cont")) %>% head() %>% cor(., use = "pairwise.complete.obs")

num_labels <- matrix(sprintf("%.2f", cor.matrix), 
                     nrow = nrow(cor.matrix), 
                     dimnames = dimnames(cor.matrix))

# diverging palette centered at 0
pheatmap(
  cor.matrix,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(200),
  breaks = seq(-1, 1, length.out = 201),   # ensures 0 is white/neutral
  display_numbers = num_labels,
  number_color = "black",
  fontsize_number = 8,
  main = "Correlation matrix (CSF, Brain, Plasma)",
  cluster_rows = F, cluster_cols = F
)

#OPCML
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_15622_13")

#CADM2
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_16907_3")

#NPTN
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_7194_36")

#GAGE2A
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_18268_5")

## GAGE2A

### Temp export forest plots
## NPTN
#CSF PC1
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_7194_36")
#CSF PC2
make.forest.plots(input.obj=raw.comp.2.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_7194_36")
#Plasma PC1
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("EDTAPlasma", test.datasets, value=T), 
                  select.prot="seq_7194_36")

## COL6A3
#CSF PC1
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_11196_31")
#CSF PC1
make.forest.plots(input.obj=raw.comp.2.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_11196_31")
#CSF PC1
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("EDTAPlasma", test.datasets, value=T), 
                  select.prot="seq_11196_31")

## TSG101
#CSF PC1
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_13044_5")
#CSF PC2
make.forest.plots(input.obj=raw.comp.2.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_13044_5")
#CSF PC1
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("EDTAPlasma", test.datasets, value=T), 
                  select.prot="seq_13044_5")


#CSF PC1
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_18268_5")
#CSF PC2
make.forest.plots(input.obj=raw.comp.2.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_18268_5")
#CSF PC1
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("EDTAPlasma", test.datasets, value=T), 
                  select.prot="seq_18268_5")



## OPCML
#CSF PC1
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_15622_13")
#CSF PC2
make.forest.plots(input.obj=raw.comp.2.lm.results, 
                  data.selection=grep("CSF", test.datasets, value=T), 
                  select.prot="seq_15622_13")
#Plasma PC1
make.forest.plots(input.obj=raw.comp.1.lm.results, 
                  data.selection=grep("EDTAPlasma", test.datasets, value=T), 
                  select.prot="seq_15622_13")


### Temp export trends
## NPTN
#CSF PC1
make.regression.plots(select.pc="Comp.1",
                      data.selection=grep("CSF", test.datasets, value=T),
                      select.prot="seq_7194_36")
#CSF PC1
make.regression.plots(select.pc="Comp.2",
                      data.selection=grep("CSF", test.datasets, value=T),
                      select.prot="seq_7194_36")
#Plasma PC1
make.regression.plots(select.pc="Comp.1",
                      data.selection=grep("EDTAPlasma", test.datasets, value=T),
                      select.prot="seq_7194_36")

## COL6A3
#CSF PC1
make.regression.plots(select.pc="Comp.1",
                      data.selection=grep("CSF", test.datasets, value=T),
                      select.prot="seq_11196_31")
#CSF PC1
make.regression.plots(select.pc="Comp.2",
                      data.selection=grep("CSF", test.datasets, value=T),
                      select.prot="seq_11196_31")
#Plasma PC1
make.regression.plots(select.pc="Comp.1",
                      data.selection=grep("EDTAPlasma", test.datasets, value=T),
                      select.prot="seq_11196_31")

## TSG101
#CSF PC1
make.regression.plots(select.pc="Comp.1",
                      data.selection=grep("CSF", test.datasets, value=T),
                      select.prot="seq_13044_5")
#CSF PC2
make.regression.plots(select.pc="Comp.2",
                      data.selection=grep("CSF", test.datasets, value=T),
                      select.prot="seq_13044_5")
#Plasma PC1
make.regression.plots(select.pc="Comp.1",
                      data.selection=grep("EDTAPlasma", test.datasets, value=T),
                      select.prot="seq_13044_5")



#### 6. Compare paired samples 

dup.ids <- person.mapping$person_id[duplicated(person.mapping$person_id)] %>% unique()
length(dup.ids) # 5512

check.dups <- metadata.gnpc.v1[metadata.gnpc.v1$sample_id %in% person.mapping$sample_id[person.mapping$person_id%in%dup.ids],] %>% 
  left_join(person.mapping[,c("sample_id", "person_id")],., by="sample_id")
table(check.dups$contributor_code, check.dups$sample_matrix) ## Solo hay dups CSF vs plasma y casi todos son nuestros

paired.samples <- check.dups %>% filter(visit==1 & sample_matrix%in%c("EDTA Plasma", "CSF")) %>% filter(person_id %in% person_id[duplicated(person_id)])
table(paired.samples$person_id) %>% quantile() # son todos duplicados

paired.samples <- paired.samples %>% 
  select(person_id, contributor_code, sample_matrix, sample_id, Dataset) %>% 
  arrange(person_id, contributor_code, sample_matrix)

### MERGE PC data:
paired.samples <- left_join(paired.samples, 
                            pcs.datasets[unique(paired.samples$Dataset)] %>% 
                              do.call('rbind',.) %>% as.data.frame() %>% 
                              tibble::rownames_to_column(var = "sample_id") %>% 
                              dplyr::select(1:6), 
                            by="sample_id")


paired.samples.wide <- paired.samples %>%
  pivot_wider(
    names_from = sample_matrix,
    values_from = -c(person_id, contributor_code),
    id_cols = c(person_id, contributor_code)
  ) %>%
  janitor::clean_names()  # replaces spaces with underscores, lowercase names, etc.


corr.matrix <- cor(paired.samples.wide[,grep("msv|comp", names(paired.samples.wide))], method = "pearson")
corr.matrix<-corr.matrix[grep("csf|msv", rownames(corr.matrix)),grep("plasma|msv", colnames(corr.matrix))]

pheatmap(corr.matrix, cluster_rows = F, cluster_cols = F, show_rownames=T, 
         show_colnames=T, display_numbers = T,breaks=seq(-1, 1, length.out=101), 
         fontsize_number = 12)


## 6. Make sorted pheatmaps

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
raw.msv.lm.results <- loadRData("Raw_lm_msv_results_20250804.RData")  # Needed to sort by msv

## Run sorted pheatmaps
lapply(test.datasets, function(x){
  
  pc.data <- pcs.datasets[[x]]
  temp.table <- preprocessed.gnpc[[x]][,-c(1:4)] 
  msv.assocs <- raw.msv.lm.results[[x]]
  
  print(paste0(x," Running sorted pheatmap; ", nrow(temp.table), " samples, ", ncol(temp.table), " proteins"))
  
  a<-pheatmap(temp.table[order(pc.data[,"msv"]),order(msv.assocs$zscore)], cluster_rows = F, cluster_cols = F, show_rownames=F, show_colnames=F)
  ggsave(paste0("Pheatmap_gradient", x, ".png"),plot=a,  width=2400, height=1700, units="px", dpi=300)
})



