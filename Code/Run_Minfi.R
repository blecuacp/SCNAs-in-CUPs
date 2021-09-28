setwd("~/your-path-to-repo/")
getwd()

rm(list=ls())

library(minfi)
library("minfiData")
library(filesstrings)

## Preprocess data
## Load 442 samples from the training set: 
file=Additional.File.2_TableS1.csv

ss <- read.csv("")

## Get methylation arrays ('.idat' files) from Biolinks organized per cancer type and per sample: 


##Biolinks code

library(TCGAbiolinks)
#Download idat files from cancer type and barcodes needed

projects <- paste("TCGA",ss$Cancer,sep = "-")
barcodes <- substr(ss$Sample_Name,1,12)

for(can in unique(projects)){
  print(can)

  Loc <- "GDCdata"
  if (file.exists(Loc)) {system(paste0("rm -r ", Loc))}

  match.file.cases.all <- NULL

  query <- tryCatch(GDCquery(project = can,
                    data.category = "Raw microarray data",
                    data.type = "Raw intensities", 
                    experimental.strategy = "Methylation array", 
                    legacy = TRUE,
                    file.type = ".idat",
                    platform = "Illumina Human Methylation 450",barcode = unique(barcodes), sample.type = "Primary solid Tumor"),
           error = function(e) query=NULL)
  
View(getResults(query))

  if(!is.null(query)){                
  match.file.cases <- getResults(query,cols=c("cases","file_name"))
  match.file.cases$project <- can
  match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
  cdir <- unique(projects[grepl(can, projects)])
  if ( ! file.exists(ccdir) ) {
    dir.create(cdir);
  }
  tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20), error = function(e) print("no files from api; method=client fails"))
  
  }
}
# This will create a map between idat file name, cases (barcode) and project
if(!is.null(query)){
  readr::write_tsv(match.file.cases.all, path =  "idat_filename_case.txt")
# code to move all files to local folder
  getwd()
for(file in dir(sprintf("GDCdata/%s/legacy/Raw_microarray_data/Raw_intensities/", can),full.names = TRUE, pattern = ".idat", recursive = T)){
  print(file)
  TCGAbiolinks:::move(from = file,to = paste(cdir,basename(file),sep = "/"))
  }

for(subdir in unique(match.file.cases.all$cases)){
  
  if ( ! dir.exists(paste(unique(match.file.cases.all$project),subdir,sep = "/"))) {
    dir.create(paste(unique(match.file.cases.all$project),subdir,sep = "/"));
  }
  
  arrays <- match.file.cases.all$file_name[which(match.file.cases.all$cases%in%subdir)]
  getwd()
  file.move(paste(unique(match.file.cases.all$project),arrays,sep = "/"), paste(unique(match.file.cases.all$project),subdir,sep = "/"),overwrite = T)
    }
  }
}


## Dump all '.idat' files all together under directory "TrainingSet_Arrays":

cdir="TrainingSet_Arrays"
if(!dir.exists(cdir))dir.create(cdir)

files <- list.files(pattern=".idat$",all.files=T,full.names=F,recursive=T)


directorio <- paste(cdir, ss$filenames,sep = "/")
rgSetn <- read.metharray(directorio,force=T)  

## quantile normalize arrays:

out <- tryCatch(genomestudion <- preprocessQuantile(rgSetn),error=function(e) {print("file-failed")})
pdf("Normalization.pdf")
densityPlot(rgSetn, sampGroups=ss$file,main="Raw", legend=FALSE)
densityPlot(getBeta(out), sampGroups=ss$file,
            main="Normalized", legend=FALSE)
dev.off()

mSetSqn <- out

## calculate detection p-values

detP <- detectionP(rgSetn)

pdf("Pvalues.pdf")
barplot(colMeans(detP), las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
dev.off()



# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSqn),rownames(detP)),]
#rm 'NA' entries
keep <- rowSums(detP < 0.01) == ncol(mSetSqn) 
table(keep)
##
mSetSqFltn <- mSetSqn[keep,]

# if your data includes males and females, remove probes on the sex chromosomes
FDATAEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- !(featureNames(mSetSqFltn) %in% FDATAEPIC$Name[FDATAEPIC$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFltn <- mSetSqFltn[keep,]

# remove probes with SNPs at CpG site
mSetSqFltn <- dropLociWithSnps(mSetSqFltn)
mSetSqFltn
?dropLociWithSnps

# exclude cross reactive probes 
library(devtools)
install_github("markgene/maxprobes", force=T)

mSetSqFltn <-  maxprobes::dropXreactiveLoci(mSetSqFltn)
mSetSqFlt <- mSetSqFltn
getwd()
save(mSetSqFlt,file = "Processed.rgSet_NSCLC")
print("Done Saving")

##16 Normal Lung Cases as controls
cdir="NormalLungControls/"
directorio <- list.files("NormalLungControls/", pattern = ".idat$")
directorio <- gsub("_Grn.idat","",directorio)
directorio <- gsub("_Red.idat","",directorio)
directorio <- unique(directorio)
directorio <- paste(cdir,directorio,sep = "/")
rgSetn <- read.metharray(directorio,force=T)  
c="NSCLC"
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

mSetSqFltn <- mSetSqFlt

getwd()
save(mSetSqFlt,file = "Processed.rgSet_TrainingSet.RData")
print("Done Saving")
