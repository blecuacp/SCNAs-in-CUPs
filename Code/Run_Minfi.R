setwd("~/your-path-to-repo/")
getwd()

rm(list=ls())

library(minfi)
library("minfiData")
library(filesstrings)

## Preprocess data
## Load 442 samples from the training set (or Validation instead): 
file="Data/Additional.File.2_TableS1.csv"

ss <- read.csv(file)

c="TrainingSet_Arrays"#ValidationSet_Arrays
directorio <- list.files(sprintf("%s",c), pattern = ".idat$", full.names = T)
directorio <- gsub("_Grn.idat","",directorio)
directorio <- gsub("_Red.idat","",directorio)
directorio <- unique(directorio)
#For Controls run:
c="WholeBlood_Controls"
directorio <- list.files(sprintf("%s",c), pattern = ".idat.gz$", full.names = T)
directorio <- gsub("_Grn.idat.gz","",directorio)
directorio <- gsub("_Red.idat.gz","",directorio)
directorio <- unique(directorio)


rgSetn <- read.metharray(directorio,force=T)  

out <- tryCatch(genomestudion <- preprocessQuantile(rgSetn),error=function(e) {print(sprintf("%s-failed",c))})
getwd()

pdf(sprintf("%s/Normalization.pdf",c))
densityPlot(rgSetn,main="Raw", legend=FALSE)
densityPlot(getBeta(out),
            main="Normalized", legend=FALSE)
dev.off()

mSetSqn <- out#genomestudion
# calculate p-values
detP <- detectionP(rgSetn)

pdf(sprintf("%s/Pvalues.pdf",c))
barplot(colMeans(detP), las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
dev.off()


## Ensure probes are in the same order in the mSetSqn and detP objects
detP <- detP[match(featureNames(mSetSqn),rownames(detP)),]

## Remove rows with at elast one 'NA' entry
keep <- rowSums(detP < 0.01) == ncol(mSetSqn) 
table(keep)
mSetSqFltn <- mSetSqn[keep,]

## If your data includes males and females, remove probes on the sex chromosomes
FDATA450 = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- !(featureNames(mSetSqFltn) %in% FDATA450$Name[FDATA450$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFltn <- mSetSqFltn[keep,]

## Remove probes with SNPs at CpG site
mSetSqFltn <- dropLociWithSnps(mSetSqFltn)
mSetSqFltn


## Exclude cross reactive probes 

mSetSqFltn <-  maxprobes::dropXreactiveLoci(mSetSqFltn)

Query <- mSetSqFltn
####For Blood controls just do: 
Controls <- mSetSqFltn
####

getwd()


save(Query,file = sprintf("%s/Processed.rgSet_%s.RData",c,c))
#For blood controls: 
save(Controls,file = sprintf("%s/Processed.rgSet_%s.RData",c,c))
print("Done Saving")
