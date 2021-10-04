## Code to download .idat files from TCGA Firehose.
## This has to be done independently for the Training and Validation Sets.
## Illustrated is the code with the TrainingSet as an example. Change input files accordingly for each case.

setwd("~/your-path-to-repo/")
getwd()

rm(list=ls())

library(minfi)
library("minfiData")
library(filesstrings)

## Load 442 samples from the training set: 
file="Additional.File.2_TableS1.csv" 

ss <- read.csv(paste("Data",file,sep = "/"))

## Get methylation arrays ('.idat' files) from Biolinks organized per cancer type and per sample within the TrainingSet_Arrays directory: 

cdir="TrainingSet_Arrays"#"ValidationSet_Arrays" for the validation ochort
if(!dir.exists(cdir))dir.create(cdir)

setwd("TrainingSet_Arrays")

##Biolinks code

library(TCGAbiolinks)
#To download idat files cancer type and barcodes needed

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
  if ( ! file.exists(cdir) ) {
    dir.create(cdir);
  }
  tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20), error = function(e) print("no files from api; method=client fails"))
  
}

# This will create a map between idat file name, cases (barcode) and project
if(!is.null(query)){
  readr::write_tsv(match.file.cases.all, path =  "idat_filename_case.txt")
# code to move all files to local folder (cancer type)
  getwd()
for(file in dir(sprintf("GDCdata/%s/legacy/Raw_microarray_data/Raw_intensities/", can),full.names = TRUE, pattern = ".idat", recursive = T)){
  print(file)
  TCGAbiolinks:::move(from = file,to = paste(cdir,basename(file),sep = "/"))
  }

# Create subdirin cancer type folder, for every barcode
  
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
setwd("../")


