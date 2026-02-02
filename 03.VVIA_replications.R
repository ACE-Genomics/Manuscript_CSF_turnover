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
library(gridExtra)

## Import ACE results
ace.results <- read.xlsx("data/Proteogenomic_signature.xlsx") %>% 
  rename(protID=AptName) %>% 
  select(protID, EntrezGeneSymbol, UniProt, Target, TargetFullName, matches("Estimate"), matches("t_value"), matches("pval"), matches("Padj_FDR")) %>% 
  select(-matches("PC"), -pval, -pval_sex, -pval_age)

names(ace.results) <- names(ace.results) %>% gsub("Estimate", "BETA", .)%>% gsub("t_value", "zscore", .) %>% gsub("pval", "P", .) %>% gsub("Padj_FDR", "FDR",.)


## Import cruchaga results:
input.dir <- "data/VVIA_replications/"
dir(input.dir) %>% dput()

variant.file.corresp <- data.frame(filename=c("forPablo_7variants_csfEURprot_excludeFACE_replication_withFeature_A1allele.csv", 
                                              "forPablo_7variants_plasmaEURmetabolomics_replication_withFeature_A1allele.csv", 
                                              "forPablo_7variants_plasmaEURproteomics_replication_withFeature_A1allele.csv", 
                                              "hg38_brain_glm_pQTL_Soma1.3k.csv", 
                                              "MQTL_brain_7SNPs.csv", 
                                              "MQTL_CSF_nonFACE_7SNPs.csv"),
                                   Tag=c("Somascan_CSF", "Metabolon_Plasma", "Somascan_Plasma","Somascan_Brain", "Metabolon_Brain", "Metabolon_CSF"))


## Read and process
variant.summary <- lapply(variant.file.corresp$filename, function(x){
  temp.table <- read.csv(paste(input.dir,x,sep="/"),h=T)
  names(temp.table)[grep("^chr", temp.table[1,])] <- "snpID" # Normalizar nomenclatura SNP
  if("LOG10_P"%in%names(temp.table)){temp.table$P <- 10^(-temp.table$LOG10_P)}  # Crear columna P when missing
  if("P"%in%names(temp.table)){temp.table$LOG10_P <- -log10(temp.table$P)}      # Crear columna LOG10P when missing
  if("Analyte"%in%names(temp.table)){names(temp.table)[names(temp.table)=="Analyte"] <- "protID"}
  if("protID"%in%names(temp.table)){temp.table$protID <- gsub("^X", "seq.", temp.table$protID)}  ## Standardize aptamer names
  if("metabID"%in%names(temp.table)){names(temp.table)[names(temp.table)=="metabID"] <- "ChemID"}
  
  temp.table
})
names(variant.summary) <- variant.file.corresp$Tag


## ADD A1 when missing (remove later)

a1.corresp <- data.frame(snpID=c("chr10:21589215:T:A", "chr11:111207001:A:G", "chr12:106083027:T:C", 
                                 "chr16:87191495:G:A", "chr22:37714443:C:T", "chr3:190902357:G:A", 
                                 "chr7:2718304:C:T"),
                         A1=c("A", "G", "C", "G", "C", "A", "T"))
variant.summary <- lapply(variant.summary, function(temp.table){
  
  if("A1"%in%names(temp.table)){temp.table}
  else{merge(a1.corresp,temp.table,by="snpID")}
  
})

lapply(variant.summary, function(temp.table){
  
  if("A1"%in%names(temp.table)){temp.table %>% select(snpID, A1) %>% unique()}
  
})

## Define ventricular vol increasing allele and flip betas

# flip.snps <- c("chr16.87191495.G.A_G", "chr7.2718304.C.T_C", "chr12.106083027.T.C_T")
define.vvias <- data.frame(snpID=c("chr10:21589215:T:A", "chr11:111207001:A:G", "chr12:106083027:T:C", 
                                   "chr16:87191495:G:A", "chr22:37714443:C:T", "chr3:190902357:G:A", 
                                   "chr7:2718304:C:T"),
                           VVIA=c("T", "G", "C", "A", "C", "G", "T"))

variant.summary <- lapply(variant.summary, function(temp.table){
  
  for(i in 1:nrow(define.vvias)){
    temp.table$BETA[temp.table$snpID==define.vvias$snpID[i]&temp.table$A1!=define.vvias$VVIA[i]] <- (-1)*temp.table$BETA[temp.table$snpID==define.vvias$snpID[i]&temp.table$A1!=define.vvias$VVIA[i]]
  }
  
  temp.table$zscore <- temp.table$BETA/temp.table$SE ## Add z score
  temp.table
  
})

### Keep CSF somascan proteins
variant.summary[["Somascan_CSF"]] <- variant.summary[["Somascan_CSF"]] %>% filter(protID%in%ace.results$protID)

### Make volcano plots: Fig 4, supp fig 8.

tested.variants <- c("chr3:190902357:G:A", "chr16:87191495:G:A", "chr7:2718304:C:T", 
                     "chr12:106083027:T:C", "chr22:37714443:C:T", "chr10:21589215:T:A", 
                     "chr11:111207001:A:G")


volcano.results<- lapply(1:length(variant.summary), function(i){
  
  temp.table <- variant.summary[[i]]
  tested.variants <- data.frame(Variant=c("chr3:190902357:G:A", "chr16:87191495:G:A", "chr7:2718304:C:T", 
                                          "chr12:106083027:T:C", "chr22:37714443:C:T", "chr10:21589215:T:A", 
                                          "chr11:111207001:A:G"),
                                Gene=c("GMNC", "C16orf95", "AMZ1/GNA12", "NUAK1", "TRIOBP", "MLLT10", "ARHGAP20/C11orf53"))
  
  volcano.list <- lapply(1:nrow(tested.variants), function(i){
    
    x <- tested.variants$Variant[i]
    
    temp.table.v1 <- temp.table[which(temp.table$snpID==x),]
    temp.table.v1$Significance <- ifelse(
      p.adjust(temp.table.v1$P, method = "fdr") < 0.05, "FDR<0.05",
      ifelse(temp.table.v1$P < 0.05, "p<0.05", "Non-significant")
    )
    temp.table.v1$Significance <- factor(
      temp.table.v1$Significance,
      levels = c("Non-significant", "p<0.05", "FDR<0.05")
    )
    
    lambda.gc <- median(qchisq(1-temp.table.v1$P,1))/qchisq(0.5,1) # guardar lambda
    
    ggplot(temp.table.v1, aes(BETA, LOG10_P, color=Significance)) + geom_point() +
      scale_color_manual(values = c("Non-significant" = "black",
                                    "p<0.05" = "blue",
                                    "FDR<0.05" = "red")) + 
      ggtitle(tested.variants$Gene[i]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
      annotate("label", x = Inf, y = Inf,
               label = paste0("?? = ", round(lambda.gc, 2)),
               hjust = 1.1, vjust = 2,
               size = 5,                   # smaller text size
               fill = alpha("white", 0.7), # semi-transparent background
               label.size = NA)     
    
  })
  
  do.call(grid.arrange, c(volcano.list, ncol = 3, top = names(variant.summary)[i]))
  
}) 
names(volcano.results) <- names(variant.summary)

## Save plots
lapply(1:length(volcano.results), function(i){
  
  ggsave(filename=paste0("VolcanoPlot_",names(volcano.results)[i], ".png"), plot=volcano.results[[i]], width=2000, height=2000, units = "px")
  
})

## Define SNPs and Genes
define.snps<-c("chr3:190902357:G:A", "chr16:87191495:G:A", "chr7:2718304:C:T", 
               "chr12:106083027:T:C", "chr22:37714443:C:T", "chr10:21589215:T:A", 
               "chr11:111207001:A:G")
define.genes <- c("GMNC", "C16orf95", "AMZ1", "NUAK1", "TRIOBP", "MLLT10", "ARHGAP20")

define.gene <- data.frame(snpID=factor(define.snps, levels=define.snps),
                          Gene=factor(define.genes, levels=define.genes),
                          orden=1:7)

## Add ACE CSF results:
variant.summary[["Somascan_ACE"]] <- ace.results %>%
  pivot_longer(
    cols = matches("^(BETA|zscore|P|FDR)\\."),   # pick only these prefixed cols
    names_to = c(".value", "Gene"),
    names_pattern = "(BETA|zscore|P|FDR)[._]?(.*)"
  ) %>% left_join(define.gene,.,by="Gene") %>% 
  select(-Gene, -orden) %>% 
  mutate(LOG10_P=-log10(P))

variant.summary <- variant.summary[c(length(variant.summary), 1:(length(variant.summary)-1))] # Reorder, ACE first

somascan.merged.tables <- lapply(variant.summary[grep("Somascan", names(variant.summary))], function(temp.table){
  
  if(!"FDR"%in%names(temp.table)){temp.table$FDR<-NA}
  merge(define.gene,temp.table,by="snpID") %>% 
    select(protID, Target, TargetFullName, UniProt, EntrezGeneSymbol, zscore, BETA, P, FDR, Gene, orden) %>% 
    arrange(orden) %>% 
    pivot_wider(
      id_cols = c(protID, EntrezGeneSymbol, UniProt, Target, TargetFullName),
      names_from = Gene,
      values_from = c(zscore, BETA, P, FDR)
    ) 
})
# names(somascan.merged.all[,grep("FDR", names(somascan.merged.all))]) <- names(variant.summary)[grep("Somascan", names(variant.summary))]

## Merge assoc stats between different tissues (pivot_wider)

somascan.merged.all <- lapply(1:length(somascan.merged.tables), function(i){
  
  temp.table <- somascan.merged.tables[[i]]
  temp.table$Tissue <- gsub("Somascan_", "", names(somascan.merged.tables)[i])
  temp.table$protID <- gsub("^X", "seq.", temp.table$protID)
  temp.table
  
}) %>% do.call('rbind', .) %>% as.data.frame() %>% 
  pivot_wider(
    id_cols = c(protID, EntrezGeneSymbol, UniProt, Target, TargetFullName),
    names_from = Tissue,
    values_from = -c(protID, EntrezGeneSymbol, UniProt, Target, TargetFullName, Tissue)
  ) 

## Check replication: FDR significant in ACE with nominal P in CSF, same direction

sapply(somascan.merged.all[,grep("FDR", names(somascan.merged.all))], function(x){sum(x<=0.05, na.rm=T)})

tissue<-"CSF"
lapply(c("FDR_GMNC_ACE", "FDR_C16orf95_ACE", "FDR_AMZ1_ACE"), function(x){
  
  gene <- gsub("FDR_", "", x) %>% gsub("_ACE", "", .)
  
  upregulated<-which(somascan.merged.all[,x]<0.05&somascan.merged.all[,paste0("BETA_",gene,"_ACE")]>0)
  downregulated<-which(somascan.merged.all[,x]<0.05&somascan.merged.all[,paste0("BETA_",gene,"_ACE")]<0)
  
  replic.upregulated <- sum(somascan.merged.all[upregulated,paste0("P_", gene, "_", tissue)]<0.05&somascan.merged.all[upregulated,paste0("BETA_", gene, "_", tissue)]>0, na.rm=T)
  replic.downregulated <- sum(somascan.merged.all[downregulated,paste0("P_", gene, "_", tissue)]<0.05&somascan.merged.all[downregulated,paste0("BETA_", gene, "_", tissue)]<0, na.rm=T)
  
  print(gene)
  print(paste0(tissue," upregulated: ", "Replicated ", replic.upregulated,"/",length(upregulated),"; ", format(100*replic.upregulated/length(upregulated), digits = 3),"%"))
  print(paste0(tissue," downregulated: ", "Replicated ", replic.downregulated,"/",length(downregulated),"; ", format(100*replic.downregulated/length(downregulated), digits = 3),"%"))
  
  
})

### z-score comparison

## Correlation heatmap all:
cor.matrix <- somascan.merged.all %>% select(matches("zscore")) %>% 
  select(matches("ACE"), matches("CSF"), matches("Brain"), matches("Plasma")) %>% 
  cor(., use="pairwise.complete.obs") 

num_labels <- matrix(sprintf("%.2f", cor.matrix), 
                     nrow = nrow(cor.matrix), 
                     dimnames = dimnames(cor.matrix))

# diverging palette centered at 0
plt4<-pheatmap(
  cor.matrix,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(200),
  breaks = seq(-1, 1, length.out = 201),   # ensures 0 is white/neutral
  display_numbers = num_labels,
  number_color = "black",
  fontsize_number = 8,
  main = "Correlation matrix (CSF, Brain, Plasma)",
  cluster_rows = F, cluster_cols = F
)

png("Heatmap_zscore_correlations_WithACE.png", width = 3250, height = 3250, res = 300)
draw(plt4)
dev.off()

## Compare the same variant across tissues
### Same variant. CSF vs brain

comparison <- c("ACE","CSF", "Brain", "Plasma")

temp.table <-lapply(comparison, function(x){
  
  merge(define.gene, variant.summary[[paste("Somascan", x, sep="_")]], by="snpID") %>% 
    mutate(Tissue = x,
           protID = gsub("^X", "seq.", protID)
    ) %>% 
    select(protID, EntrezGeneSymbol, UniProt, Target, TargetFullName, Gene, Tissue, zscore)
  
}) %>% do.call('rbind',.) %>%  pivot_wider(
  id_cols = c(protID, EntrezGeneSymbol, UniProt, Target, TargetFullName, Gene),
  names_from = Tissue,
  values_from = zscore
) 

# Compare any two columns across facets (default facet: Gene)
make_faceted_plot <- function(df, xvar, yvar, facet_var = Gene) {
  # correlation per gene
  panel_stats <- df %>%
    group_by({{ facet_var }}) %>%
    summarise(
      cor_res = list(cor.test(.data[[xvar]], .data[[yvar]], use = "complete.obs")),
      .groups = "drop"
    ) %>%
    mutate(
      r = map_dbl(cor_res, ~ .x$estimate),
      p = ifelse(map_dbl(cor_res, ~ .x$p.value)==0, "<5e-324",
                 format(map_dbl(cor_res, ~ .x$p.value), format = "e", digits = 1)),
      label = paste0("r = ", format(r, digits=2))
    )
  
  # build plot
  p <- ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(alpha = 0.8, size = 0.8, color="grey") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, color = "red") +
    facet_wrap(vars({{ facet_var }}), scales = "free_y") +
    geom_text(
      data = panel_stats,
      aes(label = label),
      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
      inherit.aes = FALSE
    ) +
    labs(
      title = paste(xvar, "vs", yvar, "z-scores"),
      x = paste(xvar, "z-score"),
      y = paste(yvar, "z-score")
    ) +
    theme_light(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "darkblue", color = "black"),
      strip.text = element_text(face = "bold")
    )
  
  # save with consistent filename
  fname <- paste0("Faceted_zscore_plot_", xvar, "vs", yvar, ".png")
  ggsave(fname, plot = p, height=1500, width=2200, units="px")
  
  return(p)
}

comparison.list <- list(c("ACE","CSF"),
                        c("ACE","Brain"),
                        c("ACE","Plasma"))

lapply(comparison.list, function(x){make_faceted_plot(temp.table, x[1], x[2])})

### EXPORT TABLE

table.for.export <- somascan.merged.all %>% 
  filter(!is.na(zscore_GMNC_ACE)) %>% 
  select(protID, EntrezGeneSymbol, UniProt, Target, TargetFullName, matches("BETA"), matches("P"), matches("FDR"), matches("zscore")) %>% 
  # select(protID, EntrezGeneSymbol, UniProt, Target, TargetFullName, matches("ACE"), matches("CSF"), matches("Brain"), matches("Plasma")) %>% 
  select(protID, EntrezGeneSymbol, UniProt, Target, TargetFullName, matches("ACE"), matches("CSF")) %>% 
  select(protID, EntrezGeneSymbol, UniProt, Target, TargetFullName, matches("GMNC"), matches("C16orf95"), matches("AMZ1"), matches("NUAK1"), 
         matches("TRIOBP"), matches("MLLT10"), matches("ARHGAP20")) %>% 
  select(where(~ !all(is.na(.)))) %>% 
  arrange(zscore_GMNC_ACE)

# write.xlsx(table.for.export,"SomaScan_WashU_replication.xlsx")


## Metabolomics 

lapply(variant.summary[grep("Metabolon", names(variant.summary))], function(x){length(unique(x$ChemID))}) 
## 1483 Plasma, 361 Brain, 440 CSF. 
lapply(variant.summary[grep("Metabolon", names(variant.summary))], function(x){unique(x[,c("CHEMICAL_NAME", "HMDB")])}) 

metabolite.list <- lapply(variant.summary[grep("Metabolon", names(variant.summary))], function(x){unique(x$ChemID)}) 

myCol <- brewer.pal(length(metabolite.list), "Pastel2")
plt<-venn.diagram(metabolite.list, category.names=grep("Metabolon", names(variant.summary), value=T), filename=NULL, fill = myCol,print.mode	=c("raw", "percent"))
grid.newpage()
grid::grid.draw(plt)


plt<-venn.diagram(metabolite.list, category.names=c("","",""), filename=NULL, fill = myCol,print.mode	=c("raw", "percent"))
grid.newpage()
grid::grid.draw(plt)


metabolon.merged.tables <- lapply(variant.summary[grep("Metabolon", names(variant.summary))], function(temp.table){
  
  temp.table$Gene<-NULL # Borro el preanotado
  merge(define.gene,temp.table,by="snpID") %>% 
    select(ChemID, CHEMICAL_NAME, HMDB, zscore, BETA, P, Gene, orden) %>% 
    arrange(orden) %>% 
    pivot_wider(
      id_cols = c(ChemID, CHEMICAL_NAME, HMDB),
      names_from = Gene,
      values_from = c(zscore, BETA, P)
    ) 
})
names(metabolon.merged.tables) <- names(variant.summary)[grep("Metabolon", names(variant.summary))]

## Merge assoc stats between different tissues (pivot_wider)
metabolon.merged.all <- lapply(1:length(metabolon.merged.tables), function(i){
  
  temp.table <- metabolon.merged.tables[[i]]
  temp.table$Tissue <- gsub("Metabolon_", "", names(metabolon.merged.tables)[i])
  temp.table
  
}) %>% do.call('rbind', .) %>% as.data.frame() %>% 
  pivot_wider(
    id_cols = c(ChemID, CHEMICAL_NAME, HMDB),
    names_from = Tissue,
    values_from = -c(ChemID, CHEMICAL_NAME, HMDB, Tissue)
  ) 


### z-score comparison (Supp. Fig. 10)

## Correlation heatmap all:
cor.matrix <- metabolon.merged.all %>% select(matches("zscore")) %>% 
  select(matches("CSF"), matches("Brain"), matches("Plasma")) %>% 
  cor(., use="pairwise.complete.obs") 

num_labels <- matrix(sprintf("%.2f", cor.matrix), 
                     nrow = nrow(cor.matrix), 
                     dimnames = dimnames(cor.matrix))

# diverging palette centered at 0
plt5 <- pheatmap(
  cor.matrix,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(200),
  breaks = seq(-1, 1, length.out = 201),   # ensures 0 is white/neutral
  display_numbers = num_labels,
  number_color = "black",
  fontsize_number = 8,
  main = "Correlation matrix (CSF, Brain, Plasma)",
  cluster_rows = F, cluster_cols = F
)


png("Heatmap_zscores_Metabolon_AllSpecies.png", width = 2600, height = 2400, res = 300)
draw(plt5)
dev.off()


