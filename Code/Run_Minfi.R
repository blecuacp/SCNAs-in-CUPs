setwd("~/your-path-to-repo/")
getwd()

rm(list=ls())

library(minfi)
library("minfiData")
library(filesstrings)

## Preprocess data
## Load 442 samples from the training set: 
file="Data/Additional.File.2_TableS1.csv"

ss <- read.csv(file)

directorio <- list.files("TrainingSet_Arrays/", pattern = ".idat$")
directorio <- gsub("_Grn.idat","",directorio)
directorio <- gsub("_Red.idat","",directorio)
directorio <- unique(directorio)
directorio <- paste(cdir,directorio,sep = "/")
rgSetn <- read.metharray(directorio,force=T)  
c="TrainigSet"
out <- tryCatch(genomestudion <- preprocessQuantile(rgSetn),error=function(e) {print(sprintf("%s-failed",c))})
getwd()

pdf("Normalization.pdf")
densityPlot(rgSetn,main="Raw", legend=FALSE)
densityPlot(getBeta(out),
            main="Normalized", legend=FALSE)
dev.off()

mSetSqn <- out#genomestudion
# calculate p-values
detP <- detectionP(rgSetn)

pdf("Pvalues.pdf")
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
FDATAEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- !(featureNames(mSetSqFltn) %in% FDATAEPIC$Name[FDATAEPIC$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFltn <- mSetSqFltn[keep,]

## Remove probes with SNPs at CpG site
mSetSqFltn <- dropLociWithSnps(mSetSqFltn)
mSetSqFltn
?dropLociWithSnps

## Exclude cross reactive probes 

mSetSqFltn <-  maxprobes::dropXreactiveLoci(mSetSqFltn)

####For Blood controls DO NOT DO THIS. Just keep and save the variable mSetSqFltn
mSetSqFltn <- mSetSqFlt
####

getwd()

#For Blood controls, save the variable mSetSqFltn
save(mSetSqFlt,file = "TrainingSet_Arrays/Processed.rgSet_TrainingSet.RData")
#For blood controls: save(mSetSqFltn,file = "WholeBlood_Controls/Processed.rgSet_TrainingSet.RData")
print("Done Saving")
