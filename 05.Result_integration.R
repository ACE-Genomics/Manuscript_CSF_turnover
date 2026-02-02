
library(openxlsx)
library(knitr)
library(data.table)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(meta)
library(grid)
library(pheatmap)

load("omics_clinical_data/Data_for_analysis_V1.RData")

pc.assoc.ace <- read.xlsx("data/Proteogenomic_signature.xlsx") 

## Integrate PC1 and PC2 GNPC associations:
pc1.assoc.gnpc <- read.table("data/Meta_Report_CSF_PC1.txt", 
                             h=T) %>% rename(AptName=apt_name)
pc2.assoc.gnpc <- read.table("data/Meta_Report_CSF_PC2.txt", 
                             h=T) %>% rename(AptName=apt_name)

names(pc1.assoc.gnpc)[-c(1:14)] <- paste(names(pc1.assoc.gnpc)[-c(1:14)], "PC1", sep=".")
names(pc2.assoc.gnpc)[-c(1:14)] <- paste(names(pc2.assoc.gnpc)[-c(1:14)], "PC2", sep=".")


rm.cols <- c("table_name", "marker", "seq_id", "seq_id_version", 
             "soma_id", "target_full_name", "target", "uni_prot", "entrez_gene_id", 
             "entrez_gene_symbol", "organism", "units", "type")

merged.pc1.assocs <- merge(pc.assoc.ace, pc1.assoc.gnpc[,!names(pc1.assoc.gnpc)%in%rm.cols],by="AptName") %>% 
  merge(.,pc2.assoc.gnpc[,!names(pc2.assoc.gnpc)%in%rm.cols],by="AptName")


## Make z-score heatmaps (supp fig. 19)

for.pheatmap <- merged.pc1.assocs 
for.pheatmap$Beta.ACE_PC1 <- for.pheatmap$Beta.Somalogic_PC.1
for.pheatmap$Beta.ACE_PC2 <- for.pheatmap$Beta.Somalogic_PC.2

pc.select <- "PC2"

lapply(c("PC1", "PC2"), function(pc.select){
  
  cor.matrix <- for.pheatmap %>% select(matches("Beta")) %>% select(matches("ACE"), matches("Cont")) %>% 
    select(matches(pc.select)) %>% select(-matches("ContQ")) %>% cor(., use="pairwise.complete.obs") 
  
  num_labels <- matrix(sprintf("%.2f", cor.matrix), 
                       nrow = nrow(cor.matrix), 
                       dimnames = dimnames(cor.matrix))
  
  # diverging palette centered at 0
  asd <- pheatmap(
    cor.matrix,
    color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(200),
    breaks = seq(-1, 1, length.out = 201),   # ensures 0 is white/neutral
    display_numbers = num_labels,
    number_color = "black",
    fontsize_number = 8,
    main = "Beta comparison",
    cluster_rows = F, cluster_cols = F
  )
  ggsave(plot = asd, paste0("Heatmap_ACE_Vs_GNPC_Beta.",pc.select,".png"), width=1600, height=1400, units="px")
  
})


## Make forest plots (supp. fig. 19)
cohorts.forest <- data.frame(Suffix=c("Somalogic", "ContJ_CSF", "ContN_CSF", "ContO_CSF", "ContT_CSF"),
                             Cohort=c("ACE", "ContJ", "ContN", "ContO", "ContT"))
for.forest <- merged.pc1.assocs
names(for.forest) <- gsub("_PC\\.", ".PC", names(for.forest))

# Define function
make.forest.plots <- function(input.obj, select.pc, select.prot){
  
  temp.meta <- data.frame(Beta=input.obj[input.obj$AptName==select.prot, paste("Beta", cohorts.forest$Suffix, select.pc, sep=".")] %>% unlist(),
                          SE=input.obj[input.obj$AptName==select.prot, paste("SE", cohorts.forest$Suffix, select.pc, sep=".")] %>% unlist(),
                          Cohort=cohorts.forest$Cohort) %>% filter(complete.cases(.))
  
  m1 <- metagen(Beta, SE, data=temp.meta, studlab=Cohort, sm="Beta")
  
  forest(m1,col.common="darkgoldenrod2",col.random="brown",col.diamond.common="darkgoldenrod2",col.diamond.random="brown",
         col.diamond.lines.random="brown",col.diamond.lines.common="darkgoldenrod2", overall = T, overall.hetstat = T, test.overall = T, 
         leftcols = c("studlab"), leftlabs = c("Cohort"), addrow = T, addrow.overall = T, colgap.forest.left = "15mm", 
         col.square = "cornflowerblue", print.tau2 = F, print.Q = T, backtransf = TRUE,scientific.pval = T, addrows.below.overall = 3,
         main="test")
  # grid.text(paste(annot.somascan$entrez_gene_symbol[annot.somascan$marker==select.prot], select.prot), .5, .9, gp=gpar(cex=1.5))
  print(paste(input.obj$EntrezGeneSymbol[input.obj$AptName==select.prot], select.prot))
}

make.forest.plots(input.obj = for.forest,
                  select.pc="PC1",
                  select.prot="seq.15622.13")

make.forest.plots(input.obj = for.forest,
                  select.pc="PC2",
                  select.prot="seq.2960.66")


### Integrate with WashU VVIA QTL reports
vvia.results <- read.xlsx("data/SomaScan_WashU_replication.xlsx")[,-c(2:5)] %>% 
  rename(AptName=protID)

names(vvia.results) <- gsub("BETA", "Estimate", names(vvia.results)) %>% gsub("^P_", "pval_", .) %>% gsub("FDR", "Padj_FDR",.) %>% gsub("zscore", "t_value",.) %>% 
  gsub("CSF$", "Cruchaga", .)

## ADD FDR col to cruchaga
cruchaga.vvias <- c("GMNC_Cruchaga", "C16orf95_Cruchaga", "AMZ1_Cruchaga", "NUAK1_Cruchaga", "TRIOBP_Cruchaga", "MLLT10_Cruchaga", "ARHGAP20_Cruchaga")

fdr.cruchaga <- lapply(cruchaga.vvias, function(x){
  p.adjust(vvia.results[,paste("pval", x, sep="_")], method = "fdr")
}) %>% do.call('cbind',.) 
colnames(fdr.cruchaga) <- paste("Padj_FDR", cruchaga.vvias, sep="_")
vvia.results <- cbind(vvia.results, fdr.cruchaga)

## Make VVIA zscore 

make.vvia.score <-  function(input.df, p.estimate="Padj_FDR", cohort, only.direction=F, only.replication=F){
  
  if(only.direction){
    sign.matrix <- input.df %>% select(matches("Estimate")) %>% select(matches(cohort)) %>% sapply(., sign) %>% rowSums() %>% return()
  }else{
    if(only.replication){
      fdr.ace <- input.df %>% select(matches("FDR")) %>% select(matches("ACE")) %>% select(matches("GMNC|C16orf95|AMZ1")) %>% sapply(., function(x){ifelse(x<0.05,1,0)})
      input.df <- input.df %>% select(matches("GMNC|C16orf95|AMZ1"))
      input.df <- input.df/fdr.ace
    }
    padj.matrix <- input.df %>% select(matches(p.estimate)) %>% select(matches(cohort)) %>% sapply(., function(x){ifelse(x<0.05,1,0)})
    sign.matrix <- input.df %>% select(matches("Estimate")) %>% select(matches(cohort)) %>% sapply(., sign)
    return(rowSums(padj.matrix*sign.matrix))
  }
  
  
}

vvia.results$VVIA_FDR_ACE <- make.vvia.score(vvia.results, p.estimate = "Padj_FDR", cohort="ACE")
vvia.results$VVIA_FDR_Cruchaga <- make.vvia.score(vvia.results, p.estimate = "Padj_FDR", cohort="Cruchaga")

vvia.results$VVIA_Nominal_ACE <- make.vvia.score(vvia.results, p.estimate = "pval", cohort="ACE")
vvia.results$VVIA_Nominal_Cruchaga <- make.vvia.score(vvia.results, p.estimate = "pval", cohort="Cruchaga")
vvia.results$VVIA_Replication_Cruchaga <- make.vvia.score(vvia.results, p.estimate = "pval", cohort="Cruchaga", only.replication=T)

vvia.results$VVIA_Direction_ACE <- make.vvia.score(vvia.results, cohort="ACE", only.direction=T)
vvia.results$VVIA_Direction_Cruchaga <- make.vvia.score(vvia.results, cohort="Cruchaga", only.direction=T)


### Make final report table (Supp. Table 9)

merged.proteogenomic.pc1 <- merge(merged.pc1.assocs[,-grep("GMNC|C16orf95|AMZ1|NUAK1|TRIOBP|MLLT10|ARHGAP20", names(merged.pc1.assocs))], 
                                  vvia.results, by="AptName")

## Select those that outperform PC1 as adjustment variables, not associated to syndromic status or disease progression.
selected.ref.proteins <- merged.proteogenomic.pc1 %>% 
  filter(
    Pbonf>0.05&ANOVA_bonf>0.05&Mean_AUC_2Y>82.96083&lambda_abeta42<10.4&lambda_ptau<14.9
  ) %>% select(AptName) %>% unlist()

merged.proteogenomic.pc1$Reference_protein <- ifelse(merged.proteogenomic.pc1$AptName%in%selected.ref.proteins, "Yes", "No")

export.table.reference.markers <- merged.proteogenomic.pc1 %>% 
  select(AptName, EntrezGeneSymbol, Reference_protein, 
         Mean_AUC_2Y, lambda_abeta42, lambda_ptau,
         HR, pval, Pbonf, ANOVA_Syndromic_Status, ANOVA_bonf, pval_sex, Pbonf_sex, pval_age, Pbonf_age,
         Beta_Random.META_PC1, Beta.Somalogic_PC.1, Beta.ContJ_CSF.PC1, Beta.ContN_CSF.PC1, Beta.ContO_CSF.PC1, Beta.ContT_CSF.PC1, 
         Beta_Random.META_PC2, Beta.Somalogic_PC.2, Beta.ContJ_CSF.PC2, Beta.ContN_CSF.PC2, Beta.ContO_CSF.PC2, Beta.ContT_CSF.PC2, 
         VVIA_FDR_ACE, VVIA_Direction_ACE, matches("t_value.*ACE"), 
         VVIA_FDR_Cruchaga,VVIA_Replication_Cruchaga, VVIA_Direction_Cruchaga, matches("t_value.*Cruchaga")) %>% 
  arrange(desc(Reference_protein),desc(Beta_Random.META_PC1))

# write.xlsx(export.table.reference.markers, "Reference_marker_Supp_Table.xlsx")

#### MAKE HEATMAP (FIGURE 4)
pc.assoc.ace <- read.xlsx("data/Proteogenomic_signature.xlsx") 

# input.data=pc.assoc.ace
# select.rows=c(pc1.top.random.effect[1:25],rev(pc1.bot.random.effect[1:25]))
# select_cols=select.cols

make.pheatmap <- function(input.data, select.rows, select.cols=NULL, cluster_rows=F, cluster_cols=F, adj.bonf=T, scale.sd=F){
  
  
  
  if(length(select.cols)==0){
    select.cols <- c("Beta.Somalogic_PC.1","Beta.Lipometrix_PC.1", "Beta_Random.META_PC1",
                     "Beta.Somalogic_PC.2", "Beta.Lipometrix_PC.2", "Beta_Random.META_PC2", 
                     "Estimate.GMNC", "Estimate.C16orf95", "Estimate.AMZ1", "Estimate.NUAK1", 
                     "Estimate.TRIOBP", "Estimate.MLLT10", "Estimate.ARHGAP20")
    
  }
  
  row.names(input.data) <- paste(input.data$EntrezGeneSymbol, input.data$AptName, sep=".")
  
  if(scale.sd){
    ref.sd <- sd(input.data[,select.cols[1]])
    input.data[,select.cols]<-sapply(input.data[,select.cols], function(x) {x*ref.sd/sd(x)})
  }
  
  for.pheatmap <- input.data[select.rows,select.cols]
  names(for.pheatmap) <- gsub("Estimate", "Beta", names(for.pheatmap))
  
  if(adj.bonf){multiple.test <- nrow(input.data)}else{multiple.test<-1}
  
  annot.pheatmap <- input.data[select.rows,] %>%
    mutate(
      Age=factor(ifelse(pval_age*multiple.test < 0.05 , "Sig.", "Non sig."), levels = c("Sig.", "Non sig.")),
      Sex=factor(ifelse(pval_sex*multiple.test < 0.05 , "Sig.", "Non sig."), levels = c("Sig.", "Non sig.")),
      Disease_prog = factor(ifelse(pval*multiple.test < 0.05 , "Sig.", "Non sig."), levels = c("Sig.", "Non sig.")),
      Disease_group = factor(ifelse(ANOVA_Syndromic_Status*multiple.test < 0.05 ,  "Sig.", "Non sig."), levels = c("Sig.", "Non sig."))
      # AUC_improve=Mean_AUC_2Y-82.26
    ) %>% 
    select(Disease_prog, Disease_group, Sex, Age)
  
  p.matrix <- input.data[select.rows,gsub("Beta", "pval", names(for.pheatmap))] %>% as.matrix()
  
  one.ast <- which(p.matrix<0.05)
  two.ast <- which(p.matrix<0.01)
  three.ast <- which(p.matrix<0.001)
  p.matrix[one.ast] <- "*"
  p.matrix[two.ast] <- "**"
  p.matrix[three.ast] <- "***"
  p.matrix[-c(one.ast,two.ast, three.ast)] <- ""
  
  annotation_colors <- list(
    Age=c("Sig." = "blue", "Non sig." = "lightgray"),
    Sex=c("Sig." = "red", "Non sig." = "lightgray"),
    Disease_prog = c("Sig." = "orange", "Non sig." = "lightgray"),
    Disease_group = c("Sig." = "purple", "Non sig." = "lightgray"))
  
  
  breaksList = c(seq(min(for.pheatmap),0,by=abs(min(for.pheatmap))/100), seq(0.01, max(for.pheatmap), by=abs(max(for.pheatmap))/100))
  
  colnames(for.pheatmap) <- gsub("Random.META_", "", colnames(for.pheatmap)) %>% paste(., "CSF", sep="_")
  
  a<-pheatmap(for.pheatmap, cluster_rows = cluster_rows, cluster_cols = cluster_cols, annotation_row = annot.pheatmap,
              display_numbers = p.matrix, 
              color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),breaks = breaksList,
              annotation_colors  = annotation_colors)
  return(a)
}

pc1.top.random.effect <- order(pc.assoc.ace$Beta_Random.META_PC1, decreasing = TRUE)[1:50]
pc1.bot.random.effect <- order(pc.assoc.ace$Beta_Random.META_PC1, decreasing = FALSE)[1:50]
pc2.top.random.effect <- order(pc.assoc.ace$Beta_Random.META_PC2, decreasing = TRUE)[1:50]
pc2.bot.random.effect <- order(pc.assoc.ace$Beta_Random.META_PC2, decreasing = FALSE)[1:50]

select.cols <- c("Beta_Random.META_PC1","Beta_Random.META_PC2", 
                 "Estimate.GMNC", "Estimate.C16orf95", "Estimate.AMZ1", "Estimate.NUAK1", 
                 "Estimate.TRIOBP", "Estimate.MLLT10", "Estimate.ARHGAP20")

a<-make.pheatmap(pc.assoc.ace, c(pc1.top.random.effect[1:25],rev(pc1.bot.random.effect[1:25])), select.cols = select.cols, scale.sd = T)
ggsave(plot = a,"Pheatmap_Comb_Top50_PC1.png", width=2000, height=3500, units="px")

a<-make.pheatmap(pc.assoc.ace, c(pc2.top.random.effect[1:25],rev(pc2.bot.random.effect[1:25])), select.cols = select.cols, scale.sd = T)
ggsave(plot = a,"Pheatmap_Comb_Top50_PC2.png", width=2000, height=3500, units="px")

## Leptin and prolactin (Supp Fig. 15)
select.genes <- c("LEP", "PRL", "TTR", "BDNF")

a <- make.pheatmap(pc.assoc.ace, which(pc.assoc.ace$EntrezGeneSymbol%in%select.genes), select.cols = select.cols, scale.sd = T)
ggsave(plot = a,"Pheatmap_LEP_PRL.png", width=3000, height=1250, units="px")


