setwd("~/Documents/CUPs/Analysis_Prem/TCGA_MSKCC/New_TCGA_Samples/Validation/")
getwd()
rm(list=ls())

load("CONUMEEvsASCAT_151_algorithm_withDiploid_AMP.AMP10_Gains_5copies.RData")
CON <- Data
load("CONUMEEvsASCAT_151_standard_with_Diploid_AMP.APM10_Gains_5copies.RData")
Ch <- Data
ss10 <- read.csv("TCGA_New_sample_IDs_With_CNVGenes_Purities_with_Diploid_5copies.csv")
dim(ss10)
ConCON=NULL
ConCh=NULL
countCoA=NULL
countChA=NULL
for(id in ss10$Sample_Name){
  
  #if(!is.null(CON$AMP[[id]]) & !is.null(Ch$AMP[[id]])){
    ConCON <- c(ConCON,CON$AMP[[id]]$Concordance,CON$AMP10[[id]]$Concordance, CON$Gain[[id]]$Concordance);
    ConCh<- c(ConCh,Ch$AMP[[id]]$Concordance);
    countCoA <- c(countCoA,CON$AMP[[id]]$NumGenes,CON$AMP10[[id]]$NumGenes, CON$Gain[[id]]$NumGenes);
    countChA <- c(countChA,Ch$AMP[[id]]$NumGenes);
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
DAT <- data.frame(CONUMEEvsASCAT_Algorithm=c(mean(ConCON)),
                  CONUMEEvsASCAT_Standard_0.3=c(mean(ConCh)),
                  CN=c(" Amp + Amp10 + Gains"))

library(reshape2)
library(ggplot2)
df = melt(DAT,id.vars = "CN")
colnames(df)[3] <- "True.Positive.Rate"
pdf(sprintf('ConcordantCalls_%s_AllGenes_Gain_General_AMp.AMP10_5copies.pdf',"AlgorythmVsStandard"),height = 5, width = 5)
print(ggplot(df, aes(CN,True.Positive.Rate , fill=variable)) + geom_bar(stat = "identity", position="dodge"))#+
     # theme(legend.position="top"))

dev.off()
getwd()

