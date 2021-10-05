## Here we store code to reproduce Fig. 1c)

setwd("~/your-path-to-repo/")
getwd()
rm(list=ls())

load("CONUMEEvsASCAT_AMP.AMP10.Gains_TPR.RData")
Algorithm<- Data

load("CONUMEEvsASCAT_Standard_TPR.RData")
Standard<- Data

file="Data/Additional.File.4_TableS3.csv"

ss <- read.csv(file)

ConAlgorithm=NULL
ConStandard=NULL

for(id in ss$Sample_Name){
  
 
    ConAlgorithm<- c(ConAlgorithm,Algorithm$AMP[[id]]$Concordance,Algorithm$AMP10[[id]]$Concordance, Algorithm$Gain[[id]]$Concordance);
    ConStandard<- c(ConStandard,Standard$AMP[[id]]$Concordance);
  
}



DAT <- data.frame(CONUMEEvsASCAT_Algorithm=c(mean(ConAlgorithm)),
                  CONUMEEvsASCAT_Standard_0.3=c(mean(ConStandard)),
                  CN=c(" Amp + Amp10 + Gains"))

library(reshape2)
library(ggplot2)
df = melt(DAT,id.vars = "CN")
colnames(df)[3] <- "True.Positive.Rate"

pdf(sprintf('ConcordantCalls_%s.pdf',"AlgorythmVsStandard"),height = 5, width = 5)
print(ggplot(df, aes(CN,True.Positive.Rate , fill=variable)) + geom_bar(stat = "identity", position="dodge"))
#+theme(legend.position="top"))
dev.off()
getwd()

