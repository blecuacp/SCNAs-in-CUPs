## Process and QC ararays, all at once. 
## This should be done independetly for both Training and Validation Sets.
## Code below with Training Set as example.

setwd("~/your-path-to-repo/")
getwd()
rm(list=ls())

library(minfi)
library("minfiData")
library(filesstrings)

## Preprocess data
cdir="WholeBlood_Controls"
files <- list.files(cdir,pattern="*.idat.gz$",all.files=T,full.names=F,recursive=F)
files <- gsub("_Grn.idat.gz","",files)
files <- gsub("_Red.idat.gz","",files)
files <- unique(files)

directory <- paste(cdir, files,sep = "/")
rgSet <- read.metharray(directory,force=T)  

## quantile normalize arrays:

out <- tryCatch(genomestudion <- preprocessQuantile(rgSet),error=function(e) {print("file-failed")})
pdf("Normalization.pdf")
densityPlot(rgSet, sampGroups=files,main="Raw", legend=FALSE)
densityPlot(getBeta(out), sampGroups=files,
            main="Normalized", legend=FALSE)
dev.off()

mSetSq <- out

## calculate detection p-values

detP <- detectionP(rgSet)

pdf("Pvalues.pdf")
barplot(colMeans(detP), las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
dev.off()



# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
#rm 'NA' entries
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)
##
mSetSqFlt <- mSetSq[keep,]

# if your data includes males and females, remove probes on the sex chromosomes
FDATA450 = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- !(featureNames(mSetSqFlt) %in% FDATA450$Name[FDATA450$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

# exclude cross reactive probes 
library(devtools)
install_github("markgene/maxprobes", force=T)

mSetSqFlt <-  maxprobes::dropXreactiveLoci(mSetSqFlt)

getwd()
save(mSetSqFlt,file = "Processed.rgSet_ControlSet.RData")
print("Done Saving")
