setwd("~/Documents/CUPs/Analysis_Prem/TCGA_MSKCC/New_TCGA_Samples/Validation/")
getwd()
rm(list=ls())

load("CONUMEEvsASCAT_151_algorithm_withDiploid_AMP.AMP10_Gains_5copies.RData")
CON <- Data
load("CONUMEEvsASCAT_151_standard_with_Diploid_AMP.APM10_Gains_5copies.RData")
Ch <- Data
load("CONUMEEvsASCAT_151_algorithm_withDiploid_FP_AMP.AMP10_Gains_5copies.RData")
CONFP <- Data
load("CONUMEEvsASCAT_151_standard_with_Diploid_FP_AMP.AMP10_5copies.RData")
ChFP <- Data
ss10 <- read.csv("TCGA_New_sample_IDs_With_CNVGenes_Purities_with_Diploid_5copies.csv")
dim(ss10)
ConCON=NULL
ConCONFP=NULL
ConCh=NULL
ConChFP=NULL
countCoA=NULL
countCoFP=NULL
countChA=NULL
countChFP=NULL
for(id in ss10$Sample_Name){
  
  #if(!is.null(CON$AMP[[id]]) & !is.null(Ch$AMP[[id]])){
    ConCON <- c(ConCON,CON$AMP[[id]]$Concordance,CON$AMP10[[id]]$Concordance,CON$Gain[[id]]$Concordance);
    ConCh<- c(ConCh,Ch$AMP[[id]]$Concordance);
    ConCONFP <- c(ConCONFP,CONFP$AMP[[id]]$Concordance,CONFP$AMP10[[id]]$Concordance,CONFP$Gain[[id]]$Concordance);
    ConChFP<- c(ConChFP,ChFP$AMP[[id]]$Concordance);
    countCoA <- c(countCoA,CON$AMP[[id]]$NumGenes,CON$AMP10[[id]]$NumGenes, CON$Gain[[id]]$NumGenes);
    countCoFP <- c(countCoFP,CONFP$AMP[[id]]$NumGenes,CONFP$AMP10[[id]]$NumGenes, CONFP$Gain[[id]]$NumGenes);
    countChA <- c(countChA,Ch$AMP[[id]]$NumGenes)
    countChFP <- c(countChFP,ChFP$AMP[[id]]$NumGenes);
  #}
  
}


#Diploid
ConCONd=NULL
ConChd=NULL
countCod=NULL
countChd=NULL
for(id in ss10$Sample_Name){
  
  #if(!is.null(CON$AMP[[id]]) & !is.null(Ch$AMP[[id]])){
  ConCONd <- c(ConCONd,CON$Diploid[[id]]$Concordance);
  ConChd<- c(ConChd,Ch$Diploid[[id]]$Concordance);
  countCod <- c(countCod,CON$Diploid[[id]]$NumGenes);
  countChd <- c(countChd,Ch$Diploid[[id]]$NumGenes);
  #}
  
}
DAT <- data.frame(CONUMEEvsASCAT_Algorithm=c(mean(ConCONFP)),
                  CONUMEEvsASCAT_Standard_0.3=c(mean(ConChFP)),
                  CN=c("FP"))

library(reshape2)
library(ggplot2)
df = melt(DAT,id.vars = "CN")
pdf(sprintf('ConcordantCalls_%s_AllGenes_AMP.AMP10_Gains_General_FP_5copies.pdf',"AlgorythmVsStandard"),height = 5, width = 5)
print(ggplot(df, aes(CN, value, fill=variable)) + geom_bar(stat = "identity", position="dodge"))
dev.off()
getwd()

### Individual Gain, Amp and Apm10

setwd("~/Documents/CUPs/Analysis_Prem/TCGA_MSKCC/New_TCGA_Samples/Validation/")
getwd()
rm(list=ls())

load("CONUMEEvsASCAT_151_algorithm_withDiploid_AMP.AMP10_Gains_5copies.RData")
CON <- Data
load("CONUMEEvsASCAT_151_standard_with_Diploid_AMP.APM10_Gains_5copies.RData")
Ch <- Data
load("CONUMEEvsASCAT_151_algorithm_withDiploid_FP_AMP.AMP10_Gains_5copies.RData")
CONFP <- Data
load("CONUMEEvsASCAT_151_standard_with_Diploid_FP_AMP.AMP10_5copies.RData")
ChFP <- Data
ss10 <- read.csv("TCGA_New_sample_IDs_With_CNVGenes_Purities_with_Diploid_5copies.csv")
dim(ss10)
ConCON=NULL
ConCONFP=NULL
for(id in ss10$Sample_Name){
  
  #if(!is.null(CON$AMP[[id]]) & !is.null(Ch$AMP[[id]])){
  ConCON <- c(ConCON,CON$AMP[[id]]$Concordance);
  #ConCh<- c(ConCh,Ch$AMP[[id]]$Concordance);
  ConCONFP <- c(ConCONFP,CONFP$AMP[[id]]$Concordance);
  #ConChFP<- c(ConChFP,ChFP$AMP[[id]]$Concordance,ChFP$HetLoss[[id]]$Concordance,ChFP$Diploid[[id]]$Concordance);
  #countCoA <- c(countCoA,CON$AMP[[id]]$NumGenes,CON$AMP10[[id]]$NumGenes,CON$Gain[[id]]$NumGenes);
  #countCoFP <- c(countCoFP,CONFP$AMP[[id]]$NumGenes,CONFP$AMP10[[id]]$NumGenes,
                 #CONFP$Gain[[id]]$NumGenes,CONFP$HetLoss[[id]]$NumGenes,CONFP$HomDel[[id]]$NumGenes,
                 #CONFP$Diploid[[id]]$NumGenes);
  #countChFP <- c(countChFP,Ch$AMP[[id]]$NumGenes,ChFP$HetLoss[[id]]$NumGenes,ChFP$Diploid[[id]]$NumGenes);
  #}
  
}

ConCON10=NULL
ConCONFP10=NULL
for(id in ss10$Sample_Name){
  
  #if(!is.null(CON$AMP[[id]]) & !is.null(Ch$AMP[[id]])){
  ConCON10 <- c(ConCON10,CON$AMP10[[id]]$Concordance);
  #ConCh<- c(ConCh,Ch$AMP[[id]]$Concordance);
  ConCONFP10 <- c(ConCONFP10,CONFP$AMP10[[id]]$Concordance);
  #ConChFP<- c(ConChFP,ChFP$AMP[[id]]$Concordance,ChFP$HetLoss[[id]]$Concordance,ChFP$Diploid[[id]]$Concordance);
  #countCoA <- c(countCoA,CON$AMP[[id]]$NumGenes,CON$AMP10[[id]]$NumGenes,CON$Gain[[id]]$NumGenes);
  #countCoFP <- c(countCoFP,CONFP$AMP[[id]]$NumGenes,CONFP$AMP10[[id]]$NumGenes,
  #CONFP$Gain[[id]]$NumGenes,CONFP$HetLoss[[id]]$NumGenes,CONFP$HomDel[[id]]$NumGenes,
  #CONFP$Diploid[[id]]$NumGenes);
  #countChFP <- c(countChFP,Ch$AMP[[id]]$NumGenes,ChFP$HetLoss[[id]]$NumGenes,ChFP$Diploid[[id]]$NumGenes);
  #}
  
}

ConCONg=NULL
ConCONFPg=NULL
for(id in ss10$Sample_Name){
  
  #if(!is.null(CON$AMP[[id]]) & !is.null(Ch$AMP[[id]])){
  ConCONg <- c(ConCONg,CON$Gain[[id]]$Concordance);
  #ConCh<- c(ConCh,Ch$AMP[[id]]$Concordance);
  ConCONFPg <- c(ConCONFPg,CONFP$Gain[[id]]$Concordance);
  #ConChFP<- c(ConChFP,ChFP$AMP[[id]]$Concordance,ChFP$HetLoss[[id]]$Concordance,ChFP$Diploid[[id]]$Concordance);
  #countCoA <- c(countCoA,CON$AMP[[id]]$NumGenes,CON$AMP10[[id]]$NumGenes,CON$Gain[[id]]$NumGenes);
  #countCoFP <- c(countCoFP,CONFP$AMP[[id]]$NumGenes,CONFP$AMP10[[id]]$NumGenes,
  #CONFP$Gain[[id]]$NumGenes,CONFP$HetLoss[[id]]$NumGenes,CONFP$HomDel[[id]]$NumGenes,
  #CONFP$Diploid[[id]]$NumGenes);
  #countChFP <- c(countChFP,Ch$AMP[[id]]$NumGenes,ChFP$HetLoss[[id]]$NumGenes,ChFP$Diploid[[id]]$NumGenes);
  #}
  
}

DAT <- data.frame(TRUE_POS=c(mean(ConCON), mean(ConCON10), mean(ConCONg)),
                  FALSE_POS=c(mean(ConCONFP),mean(ConCONFP10), mean(ConCONFPg)),
                  CN=c("Amp","Amp10","Gain"))

library(reshape2)
library(ggplot2)
df = melt(DAT,id.vars = "CN")
pdf(sprintf('ConcordantCalls_%s_AllGenes_AMP.AMP10_General_WithDiploid_FP_Amp_Amp10_Gains_5copies.pdf',"AlgorythmVsStandard"),height = 5, width = 5)
print(ggplot(df, aes(CN, value, fill=variable)) + geom_bar(stat = "identity", position="dodge"))+ theme(legend.position="top")
dev.off()
getwd()


