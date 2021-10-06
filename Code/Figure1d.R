## Here we store code to reproduce Fig. 1d)

setwd("~/your-path-to-repo/")
getwd()
rm(list=ls())


load("ValidationSet_Arrays/CONUMEEvsASCAT_AMP.AMP10.Gains_TPR.RData")
Algorithm<- Data

load("ValidationSet_Arrays/CONUMEEvsASCAT_Standard_TPR.RData")
Standard<- Data

load("ValidationSet_Arrays/CONUMEEvsASCAT_AMP.AMP10.Gains_FPR.RData")
AlgorithmFP <- Data

load("ValidationSet_Arrays/CONUMEEvsASCAT_Standard_FPR.RData")
StandardFP <- Data

file="Data/Additional.File.4_TableS3.csv"

ss <- read.csv(file)


ConAlgorithm=NULL
ConStandard=NULL
ConAlgorithmFP=NULL
ConStandardFP=NULL

for(id in ss$Sample_Name){
  
    ConAlgorithm <- c(ConAlgorithm,Algorithm$AMP[[id]]$Concordance,Algorithm$AMP10[[id]]$Concordance,Algorithm$Gain[[id]]$Concordance);
    ConStandard<- c(ConStandard,Standard$AMP[[id]]$Concordance);
    ConAlgorithmFP <- c(ConAlgorithmFP,AlgorithmFP$AMP[[id]]$Concordance,AlgorithmFP$AMP10[[id]]$Concordance,AlgorithmFP$Gain[[id]]$Concordance);
    ConStandardFP<- c(ConStandardFP,StandardFP$AMP[[id]]$Concordance)
  
}

#Figure not in the paper: FP rate Algorithm vs Standard

DAT <- data.frame(CONUMEEvsASCAT_Algorithm=c(mean(ConAlgorithmFP)),
                  CONUMEEvsASCAT_Standard_0.3=c(mean(ConStandardFP)),
                  CN=c("FP"))

library(reshape2)
library(ggplot2)
df = melt(DAT,id.vars = "CN")
pdf(sprintf('FP_%s.pdf',"AlgorythmVsStandard"),height = 5, width = 5)
print(ggplot(df, aes(CN, value, fill=variable)) + geom_bar(stat = "identity", position="dodge"))
dev.off()


### Individual Gain, Amp and Apm10
## Code for Fig. 1d)

ConAlgorithm=NULL
ConAlgorithmFP=NULL
for(id in ss$Sample_Name){
  
  
  ConAlgorithm <- c(ConAlgorithm,Algorithm$AMP[[id]]$Concordance);
  ConCONFP <- c(ConAlgorithmFP,AlgorithmFP$AMP[[id]]$Concordance);
  
}

ConAlgorithm10=NULL
ConAlgorithmFP10=NULL
for(id in ss10$Sample_Name){

  ConAlgorithm10 <- c(ConAlgorithm10,Algorithm$AMP10[[id]]$Concordance);
  ConAlgorithmFP10 <- c(ConAlgorithmFP10,AlgorithmFP$AMP10[[id]]$Concordance);
  
}

ConAlgorithmg=NULL
ConAlgorithmFPg=NULL
for(id in ss$Sample_Name){
  
  ConAlgorithmg <- c(ConAlgorithmg,Algorithm$Gain[[id]]$Concordance);
  ConAlgorithmFPg <- c(ConAlgorithmFPg,AlgorithmFP$Gain[[id]]$Concordance);
  
}

DAT <- data.frame(TRUE_POS=c(mean(ConAlgorithm), mean(ConAlgorithm10), mean(ConAlgorithmg)),
                  FALSE_POS=c(mean(ConAlgorithmFP),mean(ConAlgorithmFP10), mean(ConAlgorithmFPg)),
                  CN=c("Amp","Amp10","Gain"))

library(reshape2)
library(ggplot2)
df = melt(DAT,id.vars = "CN")
pdf(sprintf('ConcordantCalls_FP_%s.pdf',"Algorythm"),height = 5, width = 5)
print(ggplot(df, aes(CN, value, fill=variable)) + geom_bar(stat = "identity", position="dodge"))+ theme(legend.position="top")
dev.off()
getwd()


