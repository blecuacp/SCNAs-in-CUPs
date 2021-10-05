setwd("set your path/")
getwd()
rm(list=ls()) 


library(dplyr)
library(GenomicRanges)
library(regioneR)

## Load SampleSheet

file="Data/Additional.File.2_TableS1.csv"

ss <- read.csv(file)

## Load gene list

file="Data/CancerGenes.csv"
bed <- read.csv(file, stringsAsFactors = F)
dim(bed)

bed <- bed[!bed$chr%in%c("chrX","chrY"),]
dim(bed)
bed <- bed[!duplicated(bed$name),]
dim(bed)

Genes <- bed

##Load training set samples CONUMEE

file = "TrainingSet_Arrays/Log2_Ratios_TrainingSet.RData"
load(file)

##More than 10 copies
files <- list.files("path=path-to-conumee-segment-Files/",pattern="05.txt$", full.names = T, recursive = T)

## ASCAT_AMP10s: Amplifications >=10 copies


Data=list()
IDList <- NULL


for(id in files){
  
  seg <- read.delim(id, header=T,sep="\t")
  segb <- data.frame(chr=seg$chrom, start=seg$loc.start, end=seg$loc.end, log2r= seg$seg.mean)
  seggr <- toGRanges(segb)
  ID=gsub(".txt","",gsub("Segments_","",basename(id)))
  k <- which(ss$Sample_Name%in%ID)
  dd <- toGRanges(Genes)
  genes <- strsplit(as.character(ss$ASCAT_AMP10[k]),split = ";")[[1]]
  if(isEmpty(genes))next;
  int <- suppressWarnings(findOverlaps(seggr,dd[dd$name%in%genes]))
  seggr.matched <- seggr[queryHits(int)];
  mcols(seggr.matched) <- cbind.data.frame(
    mcols(seggr.matched),
    mcols(dd[dd$name%in%genes][subjectHits(int)]));
  int <- seggr.matched
  List=NULL

  if(length(int) >=1){
    IDList <- c(IDList,ID)
    for(gene in int$name){
      res <- int[int$name%in%gene,]
      List=c(List,int$log2r)
    }
    Data[[ID]]=list(y=mean(List), intercept=mean(Log2[[ID]]$log2ratio, na.rm=T),var=ss$'Purity_Impute_RFPurify(Absolute)'[ss$Sample_Name%in%ID]*sd(Log2[[ID]]$log2ratio,na.rm=T))
  }
}

length(Data)
##Fit y~x
X=NULL
for(id in IDList){
  X<-c(X,Data[[id]]$y)}
Int=NULL
for(id in IDList){
  Int<-c(Int,Data[[id]]$intercept)}

Var=NULL
for(id in IDList){
  Var<-c(Var,Data[[id]]$var)}
##Fit y~x

x <- X-Int
y <- Var

fit <-lm(Var~0+x)
summary(fit)
coeff <- fit$coefficients
K=1/coeff
getwd()
write.table(K,"Amp10_Threshold.txt", quote = F)

x <- X-Int
pdf("Amp10_fit.pdf")
plot(Var~x, ylim=c(0,1), sub=sprintf("Coeff=%.3f",coeff))
fit <-lm(Var~0+x)
abline(fit,col="red")
dev.off()
