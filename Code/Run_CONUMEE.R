setwd("set your path/")
getwd()
rm(list=ls())

library(conumee)
library(minfi)
library("minfiData")
library(regioneR)
library(filesstrings)
source("Data/Custom_Functions.R")

## Preprocess data
## Load SampleSheet Training or Validation Sets:

file="Data/Additional.File.2_TableS1.csv"

ss <- read.csv(file)

## Upload Copy Number Polymorphism file from Broad Inst. for their exclusion

file <- "Data/CNV_Germline_GSITIC2_BROAD_SNP6.merged.151117.hg19.CNV.txt"
cnvbroad <- read.delim(file,header=T,sep="\t")

## Add 'chr' tag in front of chr numbers

cnvbroad <- data.frame(chr=paste("chr",cnvbroad$Chromosome,sep = ""), start=cnvbroad$Start, end=cnvbroad$End)

## make Granges object:

cnvbroadgr <- toGRanges(cnvbroad)
Exclude <- cnvbroadgr

#Remove chr23 (X) and chr24 (Y)

Exclude <- keepSeqlevels(Exclude, paste("chr",c(1:22), sep=""), pruning.mode="coarse")
seqnames(Exclude)

## Bed file CancerGenes most amplified/most deleted in cancer autosomal chromosomes

file="Data/CancerGenes.csv"
bed <- read.csv(file, stringsAsFactors = F)
dim(bed)

bed <- bed[!bed$chr%in%c("chrX","chrY"),]
dim(bed)
bed <- bed[!duplicated(bed$name),]
dim(bed)


#Create anno object
dd <- bed
dd <- toGRanges(dd)
detail_region <- dd
seqlevels(detail_region)
anno <- CNV.create_anno(array_type = "450k",chrXY = F, 
                        exclude_regions = Exclude, detail_regions = detail_region)

dim(ss)
cancers <- unique(ss$Cancer)

c="TrainingSet_Arrays"#or ValidationSet_Arrays
##Load query samples: Training set/Validation Set
file = sprintf("%s/Processed.rgSet_%s.RData",c,c)
load(file)
my.data <- CNV.load(Query)

##Create CONUMEE/ directory:

cdir=paste(c,"CONUMEE",sep = "/")
if(!dir.exists(cdir))create_dir(cdir)
getwd()

setwd(cdir)

Log2=list()
for(can in cancers){
  print(can);
  print(which(cancers%in%can))
  
  if(!dir.exists(can))create_dir(can)
  getwd()

## Normals
## Whole Blood

## To generate the file "Processed.rgSet_BLOOD_96WB.RData" below: go to GEO repository refereed to in the paper, retrieve the '.idat' files with COde/GEO.Rand run 'Code/Run_Minfi.R' code for the 96 whole blood samples.

file="../../WholeBlood_Controls/Processed.rgSet_WholeBlood_Controls.RData"
load(file)
controls <- CNV.load(Controls)
print("control loaded")


## make sure query (my.data) and control (controls) objects have same cpg content:

my.data@intensity <- my.data@intensity[rownames(my.data@intensity)%in%rownames(controls@intensity),]
controls@intensity <-controls@intensity[rownames(controls@intensity)%in%rownames(my.data@intensity),]

for(id in ss$Sample_Name[ss$Cancer%in%can]){
  print(sprintf("muestra %i de %i %s",which(ss$Sample_Name%in%id),length(ss$Sample_Name),unique(can)));

  cdir <- paste(can,id,sep = "/")
  getwd()

  ##Make sure anno object has same #probes as control/query
  
  mset <- Query[rownames(Query)%in%rownames(my.data@intensity)]
  Mset <- mapToGenome(mset)
  anno@probes <- subsetByOverlaps(anno@probes, granges(Mset))
  
  getwd()
  if(!file.exists(cdir))dir.create(cdir)

## Run Conumee
  sample <- gsub("\\./","",ss$filenames[ss$Sample_Name%in%id])
  fit  <- CNV.fit(my.data[which(colnames(my.data@intensity)%in%sample)], controls, anno)
  Log2[[id]] <- list(log2ratio=fit@fit$ratio, Purity=ss$Purity_Impute_RFPurify.Absolute.[ss$Sample_Name%in%id])
  fit2 <- CNV.segment(CNV.detail(CNV.bin(fit)))

##Visualize and write up coordinates
getwd()
pdf(file = sprintf("%s_annotated.pdf",id),height = 9, width = 18)
CNV.genomeplot2(fit2)
dev.off()

write.table(CNV.write(fit2, what = "segments"),sprintf("Segments_%s.txt",id),col.names = T, row.names = F, quote = F, sep = "\t") 
file.move(c(sprintf("%s_annotated.pdf",id),sprintf("Segments_%s.txt",id)),cdir, overwrite = T)
  }
}

save(Log2, file=sprintf("Log2_Ratios_%s.RData",c))

setwd("../../")
