setwd("/set-your-path/")
getwd()
rm(list=ls())

#CONUMEE Calibration with TCGAvsBLOOD: K_c cnostants from the text:
K_amp10= 3.674
K_amp = 1.54
K_gain= 1.118
K_hetloss=1.885
K_homdel=4.9

## Uplodad Validation set TableS3:

file="Data/Additional.File.4_TableS3.csv"
ss <- read.csv(file)
dim(ss)

##Run 'Run_Minfi.R' and 'Run_CONUMEE.R' for validation set above to generate: "Log2_Ratios_ValidationSet_Arrays.RData"

load("ValidationSet_Arrays/CONUMEE/Log2_Ratios_ValidationSet_Arrays.RData")
getwd()
IDs <- list.files(path="ValidationSet_Arrays/CONUMEE/", pattern = "^Segments", full.names = TRUE, recursive = TRUE)

Dirs <- gsub("Segments_","",gsub(".txt","",basename(IDs)))

head(Dirs)

for(id in Dirs){
  if(is.null(Log2[[id]]$log2ratio))next;
  Log2[[id]]$Gain=mean(Log2[[id]]$log2ratio, na.rm=T) + ss$Purity_Impute_Absolute[which(ss$Sample_Name%in%id)]*K_gain*sd(Log2[[id]]$log2ratio, na.rm = T)
  Log2[[id]]$Amp=mean(Log2[[id]]$log2ratio, na.rm = T) + ss$Purity_Impute_Absolute[which(ss$Sample_Name%in%id)]*K_amp*sd(Log2[[id]]$log2ratio,na.rm = T)
  Log2[[id]]$Amp10=mean(Log2[[id]]$log2ratio, na.rm = T) + ss$Purity_Impute_Absolute[which(ss$Sample_Name%in%id)]*K_amp10*sd(Log2[[id]]$log2ratio,na.rm = T)
  Log2[[id]]$HetLoss=mean(Log2[[id]]$log2ratio, na.rm = T) - ss$Purity_Impute_Absolute[which(ss$Sample_Name%in%id)]*K_hetloss*sd(Log2[[id]]$log2ratio,na.rm = T)
  Log2[[id]]$HomDel=mean(Log2[[id]]$log2ratio, na.rm = T) - ss$Purity_Impute_Absolute[which(ss$Sample_Name%in%id)]*K_homdel*sd(Log2[[id]]$log2ratio,na.rm = T)
  
}

##Threshold calls

for(dir in Dirs){
  
  print(which(Dirs%in%dir))
  cancer=ss$Cancer[ss$Sample_Name%in%dir]
  seg <- read.delim(sprintf("ValidationSet_Arrays/CONUMEE/%s/%s/Segments_%s.txt",cancer,dir,dir), header=T, sep="\t")
  dim(seg)
  if(is.null(Log2[[dir]]$log2ratio))next;
  seg$CNCall="Diploid"
  for(i in 1:nrow(seg)){
    if(seg$seg.mean[i] >= Log2[[dir]]$Gain )seg$CNCall[i] <- "Gain"
    if(seg$seg.mean[i] >= Log2[[dir]]$Amp )seg$CNCall[i] <- "Amp"
    if(seg$seg.mean[i] >= Log2[[dir]]$Amp10 )seg$CNCall[i] <- "Amp10"
    if(seg$seg.mean[i] <= Log2[[dir]]$HetLoss )seg$CNCall[i] <- "HetLoss"
    if(seg$seg.mean[i] <= Log2[[dir]]$HomDel )seg$CNCall[i] <- "HomDel"
  }
  write.table(seg,sprintf("ValidationSet_Arrays/CONUMEE/%s/%s/Segments_%s_withCalls.txt",cancer,dir,dir), col.names = T, row.names = F, quote = F,sep="\t")
  
}

### Standard thresholding: +/-0.3

for(dir in Dirs){

  print(which(Dirs%in%dir))
  cancer=ss$Cancer[ss$Sample_Name%in%dir]
  seg <- read.delim(sprintf("ValidationSet_Arrays/CONUMEE/%s/%s/Segments_%s.txt",cancer,dir,dir), header=T, sep="\t")
  dim(seg)
  if(is.null(Log2[[dir]]$log2ratio))next;
  seg$CNCall="Diploid"

## >= 0.3 ==> Amp; < -0.3 ===> HetLoss
  for(i in 1:nrow(seg)){
    if(seg$seg.mean[i] >= 0.3)seg$CNCall[i] <- "Amp"
    if(seg$seg.mean[i] <= -0.3)seg$CNCall[i] <- "HetLoss"
  }
  write.table(seg,sprintf("ValidationSet_Arrays/CONUMEE/%s/%s/Segments_%s_withCalls_Standard.txt",cancer,dir,dir), col.names = T, row.names = F, quote = F,sep="\t")
}
