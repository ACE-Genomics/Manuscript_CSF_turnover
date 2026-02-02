#!/usr/bin/env Rscript

args=commandArgs(trailingOnly = TRUE)

library(data.table)
library(reshape2)

sumstats.dir <- args[1]
freq.variants.path <- args[2]
output <- args[3]

### STEP 1: Read Freq file.
freq.variants <- fread(freq.variants.path, h=T)
names(freq.variants)<- gsub("\\#CHROM", "CHR", names(freq.variants))

### STEP 2: Process and align individual sumstats

all.sumstats<-dir(sumstats.dir)

keep.cols <- c("ID", "OBS_CT", "BETA", "SE", "T_STAT", "P")
linear.flip <- c("BETA", "T_STAT")

result.sumstats <- list()

for(i in 1:length(all.sumstats)){
  
  sumstats.path <- paste(sumstats.dir, "/", all.sumstats[i], "/", all.sumstats[i], ".glm.linear", sep="")
  print(paste("READING SUMSTATS:",sumstats.path))
  sumstats<- fread(sumstats.path, h=T)
  flip <- which(sumstats$A1==sumstats$REF)
  
  sumstats[flip,"BETA"] <- (-1)*sumstats[flip,"BETA"]
  sumstats[flip,"T_STAT"] <- (-1)*sumstats[flip,"T_STAT"]
  
  sumstats.v1<-sumstats[, .(ID, OBS_CT, BETA, SE, T_STAT, P)]
  names(sumstats.v1)[-1] <- paste(names(sumstats.v1)[-1], all.sumstats[i], sep=".")
  
  result.sumstats[[i]] <- sumstats.v1
}


### STEP 3: Fast merge for data.tables
print("Merging sumstats")
merged.sumstats <-Reduce(function(x, y) x[y, on = "ID"], result.sumstats)
merged.sumstats <- merge(freq.variants,merged.sumstats,by="ID")
# merged.sumstats<-cbind(colsplit(merged.sumstats$ID, pattern = ":", names = c("CHR", "POS", "REF", "ALT")),merged.sumstats)

### Write output
fwrite(merged.sumstats, output, sep="\t")
