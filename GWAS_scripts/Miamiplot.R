#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(miamiplot)

# gwas.a <- "/nas/HARPONE/QTL_Analysis/mQTL/03.Retomar_trabajo_apr2024/GWAS/Lipometrix_PC.1/Lipometrix_PC.1.glm.linear"
# gwas.b <- "/nas/HARPONE/QTL_Analysis/mQTL/03.Retomar_trabajo_apr2024/GWAS/Somalogic_PC.1/Somalogic_PC.1.glm.linear"
# out.dir <- "/nas/HARPONE/QTL_Analysis/mQTL/03.Retomar_trabajo_apr2024/report_GWAS"
# out.prefix <- "Lipidomics_Proteomics_PC1"

gwas.a <- args[1]
gwas.b <- args[2]
out.dir <- args[3]
out.prefix <- args[4]

######

## Sumstats A
# Leer sumstats
for.miami.a <- fread(gwas.a, h=T)

# Adaptar sumstats:
for.miami.a <- for.miami.a[,c("CHROM", "POS", "ID","P",  "ERRCODE","A1_FREQ")]
names(for.miami.a) <- c("CHR", "BP", "SNP", "P", "ERRCODE","A1_FREQ")

for.miami.a$P <- as.numeric(for.miami.a$P)
for.miami.a <- for.miami.a[!is.na(for.miami.a$P),]
for.miami.a <- for.miami.a[which(for.miami.a$ERRCODE=="."),]
for.miami.a <- for.miami.a[,c("CHR", "BP", "SNP", "P")]
for.miami.a$Study <- "A"

## Sumstats B
# Leer sumstats
for.miami.b <- fread(gwas.b, h=T)

# Adaptar sumstats:
for.miami.b <- for.miami.b[,c("CHROM", "POS", "ID","P",  "ERRCODE","A1_FREQ")]
names(for.miami.b) <- c("CHR", "BP", "SNP", "P", "ERRCODE","A1_FREQ")

for.miami.b$P <- as.numeric(for.miami.b$P)
for.miami.b <- for.miami.b[!is.na(for.miami.b$P),]
for.miami.b <- for.miami.b[which(for.miami.b$ERRCODE=="."),]
for.miami.b <- for.miami.b[,c("CHR", "BP", "SNP", "P")]
for.miami.b$Study <- "B"

## Combine sumstats and make plot
combined.sumstats <- rbind(for.miami.a,for.miami.b)

miami.png <- paste(out.dir, "/", "Miamiplot_", out.prefix,".png", sep="")
png(miami.png, 500, 400)
ggmiami(combined.sumstats, split_by="Study", split_at = "A",  chr="CHR", pos="BP", p="P",
        upper_ylab = gsub("_", " ", gsub(".glm.linear|\\.","",basename(gwas.a))), lower_ylab = gsub("_", " ", gsub(".glm.linear|\\.","",basename(gwas.b))))

dev.off()
