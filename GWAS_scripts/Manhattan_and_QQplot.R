args=commandArgs(trailingOnly = TRUE)

library(data.table)
library(qqman)

input.data <- args[1]
out.dir <- args[2]
out.prefix <- args[3]

# Leer sumstats
for.manhattan <- fread(input.data, h=T)

# Adaptar sumstats:
for.manhattan <- for.manhattan[,c("CHROM", "POS", "ID","P",  "ERRCODE","A1_FREQ", "OBS_CT")]
names(for.manhattan) <- c("CHR", "BP", "SNP", "P", "ERRCODE","A1_FREQ", "OBS_CT")

for.manhattan$P <- as.numeric(for.manhattan$P)
for.manhattan <- for.manhattan[!is.na(for.manhattan$P)]
for.manhattan <- for.manhattan[which(for.manhattan$ERRCODE=="."),]
keep.maf05 <- which(for.manhattan$A1_FREQ>=0.05&for.manhattan$A1_FREQ<=0.95)

n.snps <- nrow(for.manhattan) # GUardar N SNPs report
n.samples <- mean(for.manhattan$OBS_CT)


# Info
print(paste("PASS SNPs:", nrow(for.manhattan)))
print(paste("PASS SNPs (MAF>=0.05):", length(keep.maf05)))

for.manhattan <- for.manhattan[,c("CHR", "BP", "SNP", "P")]

# Make Manhattan Plot
manhattan.png <- paste(out.dir, "/", "Manhattan.", out.prefix,".png", sep="")
#manhattan.title <- paste(plot.title)

png(manhattan.png,1000, 500)
manhattan(for.manhattan, col = c("azure4", "gray7"))
dev.off()

# Make QQPlot (MAF>0.05)

for.qqplot<-for.manhattan[keep.maf05,]
lambda.gc <- median(qchisq(1-for.qqplot$P,1))/qchisq(0.5,1) # guardar lambda

qqplot.png <- paste(out.dir, "/", "QQPlot_MAF0.05", out.prefix,".png", sep="")

png(qqplot.png,750, 750)
par(mar=c(5,6,4,2)+.1) # AUMENTAR VENTANA PARA QUE ENTREN LOS XLABs
qq(for.qqplot$P, cex=1.5,cex.lab=2,cex.axis=2)
legend("topleft",paste("λGC=",round(lambda.gc,digits=3),sep=""), cex=2)
dev.off()

# Summary report: 
report <- data.frame(SPECIES=out.prefix,N_SAMPLES=n.samples,N_SNPs=n.snps,lambdaGC=lambda.gc)
write.table(report, paste(out.dir, "/", "GWAS_report_", out.prefix,".txt", sep=""), row.names=F,col.names=F, quote=F, sep="\t")


