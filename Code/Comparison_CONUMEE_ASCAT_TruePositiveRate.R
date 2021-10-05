## Here we calculate the True Positive Rate (TPR) for each amplification case c={"Amp","Amp10","Gain"}, 
## on a per sample basis. We use the Validation Set only.


setwd("~/your-path-to-repo/")
getwd()
rm(list=ls())

## Using our calibration method



#CONUMEE Calibration with TCGAvsBLOOD
amp10= 3.674
amp = 1.54
gain= 1.118
hetloss=1.885
homdel=4.9


file="Data/CancerGenes.csv"
bed <- read.csv(file, stringsAsFactors = F)
dim(bed)

bed <- bed[!bed$chr%in%c("chrX","chrY"),]
dim(bed)
bed <- bed[!duplicated(bed$name),]
dim(bed)

file="Data/Additional.File.4_TableS3.csv"

ss <- read.csv(file)

##TPRs: True Positive Rates, defined as TP/(TP+FN)
##Amplification: "Amp" called by CONUMEE that are "Amp" in TCGA

files <- list.files(path = "CONUMEE_ValidationSet/",pattern = "Calls.txt$",full.names = T,recursive = T)
Data=list()

for(id in files){
  print(which(files%in%id))
  seg <- read.delim(id, header=T,sep="\t")
  segb <- data.frame(chr=seg$chrom, start=seg$loc.start, end=seg$loc.end, log2r= seg$seg.mean, CNCall=seg$CNCall)
  seggr <- toGRanges(segb)
  ID=gsub("_withCalls.txt","",gsub("Segments_","",basename(id)))
  k <- which(ss10$Sample_Name%in%ID)
  dd <- toGRanges(bed)
  genes <- strsplit(as.character(ss$ASCAT_AMP),split = ";")[[1]]
  int <- suppressWarnings(findOverlaps(seggr,dd[dd$name%in%genes]))
  seggr.matched <- seggr[queryHits(int)];
  # Add the metadata from dd
  mcols(seggr.matched) <- cbind.data.frame(
    mcols(seggr.matched),
    mcols(dd[dd$name%in%genes][subjectHits(int)]));
  #i=1
  count=0
  int <- seggr.matched[!duplicated(seggr.matched$name)]
  if(length(int) >=1){
    for(i in 1:length(int)){
      #print(i)
      res <- int[i,]
      dat <- filter(segb, chr%in%as.character(res@seqnames@values) & start <= res@ranges@start & end >= (res@ranges@start + res@ranges@width-1) & CNCall=="Amp")
      if(nrow(dat)!=0)count=count+1;
    }
  Data[["AMP"]][[ID]]=list(Count=count, NumGenes=length(unique(seggr.matched$name)),Concordance=count/length(unique(seggr.matched$name)))
  }
}

## in Data[["AMP"]][[ID]] --> for every ID: Count==TP; FN==NumGenes-TP; Concordance==TPR=TP/(TP+FN)=TP/NumGenes.
## Same holds below for Amp10 and Gain.

##Amplifications10: "Amp10" called by CONUMEE that are "Amp10" in TCGA
for(id in files){
  print(which(files%in%id))
  seg <- read.delim(id, header=T,sep="\t")
  segb <- data.frame(chr=seg$chrom, start=seg$loc.start, end=seg$loc.end, log2r= seg$seg.mean, CNCall=seg$CNCall)
  seggr <- toGRanges(segb)
  ID=gsub("_withCalls.txt","",gsub("Segments_","",basename(id)))
  k <- which(ss10$Sample_Name%in%ID)
  dd <- toGRanges(bed)
  genes <- strsplit(as.character(ss$ASCAT_AMP10[k]),split = ";")[[1]]
  int <- suppressWarnings(findOverlaps(seggr,dd[dd$name%in%genes]))
  seggr.matched <- seggr[queryHits(int)];
  # Add the metadata from dd
  mcols(seggr.matched) <- cbind.data.frame(
    mcols(seggr.matched),
    mcols(dd[dd$name%in%genes][subjectHits(int)]));
  #i=1
  count=0
  int <- seggr.matched[!duplicated(seggr.matched$name)]
  if(length(int) >=1){
    for(i in 1:length(int)){
      #print(i)
      res <- int[i,]
      dat <- filter(segb, chr%in%as.character(res@seqnames@values) & start <= res@ranges@start & end >= (res@ranges@start + res@ranges@width-1) & CNCall=="Amp10")
      if(nrow(dat)!=0)count=count+1;
    }
    Data[["AMP10"]][[ID]]=list(Count=count, NumGenes=length(unique(seggr.matched$name)),Concordance=count/length(unique(seggr.matched$name)))
  }
}

##Gains:"Gain" called by CONUMEE that are "Gain" in TCGA 
for(id in files){
  print(which(files%in%id))
  seg <- read.delim(id, header=T,sep="\t")
  segb <- data.frame(chr=seg$chrom, start=seg$loc.start, end=seg$loc.end, log2r= seg$seg.mean, CNCall=seg$CNCall)
  seggr <- toGRanges(segb)
  ID=gsub("_withCalls.txt","",gsub("Segments_","",basename(id)))
  #ID="TCGA-05-4433-01A-22D-1856-05"
  k <- which(ss10$Sample_Name%in%ID)
  dd <- toGRanges(bed)
  genes <- strsplit(as.character(ss$ASCAT_Gain[k]),split = ";")[[1]]
  int <- suppressWarnings(findOverlaps(seggr,dd[dd$name%in%genes]))
  seggr.matched <- seggr[queryHits(int)];
  # Add the metadata from dd
  mcols(seggr.matched) <- cbind.data.frame(
    mcols(seggr.matched),
    mcols(dd[dd$name%in%genes][subjectHits(int)]));
  #i=1
  count=0
  int <- seggr.matched[!duplicated(seggr.matched$name)]
  if(length(int) >=1){
    for(i in 1:length(int)){
      #print(i)
      res <- int[i,]
      dat <- filter(segb, chr%in%as.character(res@seqnames@values) & start <= res@ranges@start & end >= (res@ranges@start + res@ranges@width-1) & CNCall=="Gain")
      if(nrow(dat)!=0)count=count+1;
    }
    Data[["Gain"]][[ID]]=list(Count=count, NumGenes=length(unique(seggr.matched$name)),Concordance=count/length(unique(seggr.matched$name)))
  }
}


save(Data,file="CONUMEEvsASCAT_AMP.AMP10.Gains_TPR.RData")

## Using the standard method: 

### if log2(R) >0.3 ==> "Amp" (All gains,amps and amp10s in our algorithm) and "Amp" or "Amp10" or "Gains" in TCGA 

files <- list.files(path = "CONUMEE_ValidationSet/",pattern = "Calls_Standard.txt$",full.names = T,recursive = T)
Data=list()

for(id in files){
  print(which(files%in%id))
  seg <- read.delim(id, header=T,sep="\t")
  segb <- data.frame(chr=seg$chrom, start=seg$loc.start, end=seg$loc.end, log2r= seg$seg.mean, CNCall=seg$CNCall)
  seggr <- toGRanges(segb)
  ID=gsub("_withCalls_Standard.txt","",gsub("Segments_","",basename(id)))
  k <- which(ss10$Sample_Name%in%ID)
  dd <- toGRanges(bed)
  genes <- c(strsplit(as.character(ss$ASCAT_AMP[k]),split = ";")[[1]],
             strsplit(as.character(ss$ASCAT_AMP10[k]),split = ";")[[1]],
             strsplit(as.character(ss$ASCAT_Gain[k]),split = ";")[[1]])
  #genes <- unique(genes)
  int <- suppressWarnings(findOverlaps(seggr,dd[dd$name%in%genes]))
  seggr.matched <- seggr[queryHits(int)];
  # Add the metadata from dd
  mcols(seggr.matched) <- cbind.data.frame(
    mcols(seggr.matched),
    mcols(dd[dd$name%in%genes][subjectHits(int)]));
  #i=1
  count=0
  int <- seggr.matched[!duplicated(seggr.matched$name)]
  if(length(int) >=1){
    for(i in 1:length(int)){
      #print(i)
      res <- int[i,]
      dat <- filter(segb, chr%in%as.character(res@seqnames@values) & 
                      start <= res@ranges@start & end >= (res@ranges@start + res@ranges@width-1) & 
                      CNCall%in%c("Amp"))
      if(nrow(dat)!=0)count=count+1;
    }
    Data[["AMP"]][[ID]]=list(Count=count, NumGenes=length(unique(seggr.matched$name)),Concordance=count/length(unique(seggr.matched$name)))
  }
}


### if log2(R) <-0.3 ==> "HetLoss" (HetLosses and HomDels in our algorithm) and "HetLoss" or "Homdel" in TCGA


for(id in files){
  print(which(files%in%id))
  seg <- read.delim(id, header=T,sep="\t")
  segb <- data.frame(chr=seg$chrom, start=seg$loc.start, end=seg$loc.end, log2r= seg$seg.mean, CNCall=seg$CNCall)
  seggr <- toGRanges(segb)
  ID=gsub("_withCalls_Standard.txt","",gsub("Segments_","",basename(id)))
  k <- which(ss10$Sample_Name%in%ID)
  dd <- toGRanges(bed)
  genes <- c(strsplit(as.character(ss$ASCAT_HetLoss[k]),split = ";")[[1]],
             strsplit(as.character(ss$ASCAT_HomDel[k]),split = ";")[[1]])
             
  int <- suppressWarnings(findOverlaps(seggr,dd[dd$name%in%genes]))
  seggr.matched <- seggr[queryHits(int)];
  # Add the metadata from dd
  mcols(seggr.matched) <- cbind.data.frame(
    mcols(seggr.matched),
    mcols(dd[dd$name%in%genes][subjectHits(int)]));
  #i=1
  count=0
  int <- seggr.matched[!duplicated(seggr.matched$name)]
  if(length(int) >=1){
    for(i in 1:length(int)){
      #print(i)
      res <- int[i,]
      dat <- filter(segb, chr%in%as.character(res@seqnames@values) & 
                      start <= res@ranges@start & end >= (res@ranges@start + res@ranges@width-1) & 
                      CNCall%in%c("HetLoss"))
      if(nrow(dat)!=0)count=count+1;
    }
    Data[["HetLoss"]][[ID]]=list(Count=count, NumGenes=length(unique(seggr.matched$name)),Concordance=count/length(unique(seggr.matched$name)))
  }
}


save(Data,file="CONUMEEvsASCAT_Standard_TPR.RData")
