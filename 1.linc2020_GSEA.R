require(data.table)
require(stringi)
require(parallel)
library(data.table)
library(cmapR)
library(dplyr)
rm(list=ls())
setwd("workingDirectory/")
#function
source('function.R')
#1.metadata#########################################################################################################################
meta=fread("instinfo_beta.txt")
metadata=meta[,c("rna_plate",
  "rna_well",  
  "pert_id",
  "cmap_name",
  "pert_type",
  "pert_dose",
  "pert_dose_unit",
  "pert_time",
  "pert_time_unit",
  "cell_mfc_name",
  "sample_id")]
metadata$pert_iname=metadata$cmap_name
metadata$cell_id=metadata$cell_mfc_name
metadata$Num=metadata$sample_id


# Global setting
# Output
inDir = c('Data'='Data')
outDir = './Data/Compound_Data/'
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
metadata$pert_cdose = sprintf(
    '%s%s',
    metadata$pert_dose,
    metadata$pert_dose_unit
)
metadata$pert_cdose[grepl('^NA', metadata$pert_cdose)] = ''
metadata$pert_ctime = sprintf(
    '%s%s',
    metadata$pert_time,
    metadata$pert_time_unit
)
metadata$pert_cval = sprintf(
    '%s %s for %s in %s',
    metadata$pert_cdose,
    metadata$pert_iname,
    metadata$pert_ctime,
    metadata$cell_id
)
metadata$pert_cval = stri_trim_left(metadata$pert_cval)
metadata$rna_centre = sapply(strsplit(x = as.character(metadata$rna_plate), split = '_'), 
                               USE.NAMES = FALSE, FUN = function(i) {return(i[[1]])})
metadata= subset(metadata, pert_type %in% c('ctl_vehicle', 'ctl_untrt','trt_cp'))
sampleMetadata = metadata[,c(
    'Num',
    'rna_plate',
    'rna_centre',
    'pert_iname',
    'pert_type',
    'pert_dose',
    'pert_dose_unit',
    'pert_time',
    'pert_cdose',
    'pert_ctime',
    'pert_cval',
    'cell_id'
)]
sampleMetadata$source_dataset = rep("Data", nrow(metadata))



# Variable Cleanup
rm( metadata)
# ----- Generate Baseline-Matched Condition Metadata -----
cat('Generating Baseline-Matched Condition Metadata.\n')
  # Match Contrast-Baseline
selectionVector = c(
    'source_dataset',
    'pert_cval',
    'rna_centre',
    'pert_type',
    'cell_id',
    'pert_iname',
    'pert_ctime',
    'pert_cdose'
)
useDF = unique(sampleMetadata[, c(
    'source_dataset',
    'pert_cval',
    'rna_centre',
    'pert_type',
    'cell_id',
    'pert_iname',
    'pert_ctime',
    'pert_cdose'
)])
colnames(useDF)

#Sample Counts
sampleMetadata = as.data.table(sampleMetadata)
setkey(sampleMetadata, rna_centre, pert_type, pert_cval)
useDF$sampleCount = unlist(mclapply(1:nrow(useDF), mc.preschedule = TRUE, mc.cores = 160, mc.cleanup = TRUE, FUN = function(i) {
    tempCondition = useDF[i, ]
    return(sampleMetadata[.(tempCondition$rna_centre, tempCondition$pert_type, tempCondition$pert_cval), .N, nomatch = 0])
  }), recursive = FALSE, use.names = FALSE)
sampleMetadata = as.data.frame(sampleMetadata)
colnames( useDF)
# Filter out Sample Size == 1
useDF = subset(useDF, sampleCount > 1)
# Prepare Subsets
casesDF = subset(useDF, grepl('^trt', useDF$pert_type))
controlsDF = subset(useDF, grepl('^ctl', useDF$pert_type))
summary(casesDF)  
#function
SelectionResolver = function(conditionDF, controlDF) {
  # Output Resolving Function
  finalSampleSize = min(conditionDF$sampleCount, controlDF$sampleCount, 20)
  tempVector = c(
    'source_dataset' = conditionDF$source_dataset,
    'rna_centre' = conditionDF$rna_centre,
    'cell_id' = conditionDF$cell_id,
    
    'trt_type' = conditionDF$pert_type,
    'trt_cval' = conditionDF$pert_cval,
    'trt_iname' = conditionDF$pert_iname,
    'trt_cdose' = conditionDF$pert_cdose,
    'trt_ctime' = conditionDF$pert_ctime,
    'trt_orig_sampleCount' = conditionDF$sampleCount,
    'trt_final_sampleCount' = finalSampleSize,
    
    'ctl_type' = controlDF$pert_type,
    'ctl_cval' = controlDF$pert_cval,
    'ctl_iname' = controlDF$pert_iname,
    'ctl_cdose' = controlDF$pert_cdose,
    'ctl_ctime' = controlDF$pert_ctime,
    'ctl_orig_sampleCount' = controlDF$sampleCount,
    'ctl_final_sampleCount' = finalSampleSize
  )
  return(tempVector)
}

tempList = mclapply(1:nrow(casesDF), mc.preschedule = TRUE, mc.cores = 10, mc.cleanup = TRUE, FUN = function(currentRowIndex) {
    condition_aaa = casesDF[currentRowIndex, ]
    centre_bbb = condition_aaa$rna_centre
    type_ccc = condition_aaa$pert_type
    cell_ddd = condition_aaa$cell_id
    time_eee = condition_aaa$pert_ctime
    dose_fff = condition_aaa$pert_cdose
    
    controlSubset = subset(controlsDF, rna_centre == centre_bbb & cell_id == cell_ddd & pert_ctime == time_eee)
    
    # Compound Data Resolving
    if (type_ccc == 'trt_cp') {
      controlSubset = subset(controlSubset, pert_type %in% c('ctl_vehicle', 'ctl_untrt'))
      
      # SKIP if No Matching: 
      if (nrow(controlSubset) == 0) {
        return(NULL)
      }
      
      if ('DMSO' %in% controlSubset$pert_iname) {
     
        controlSubset = subset(controlSubset, pert_iname == 'DMSO')
        finalOutput = controlSubset[order(controlSubset$sampleCount, decreasing = TRUE), ][1, ]
        return(SelectionResolver(conditionDF = condition_aaa, controlDF = finalOutput))
      }
      
      if (any(c('PBS', 'H2O', 'UnTrt') %in% controlSubset$pert_iname)) {
        
        finalOutput = controlSubset[order(controlSubset$sampleCount, decreasing = TRUE), ][1, ]
        return(SelectionResolver(conditionDF = condition_aaa, controlDF = finalOutput))
        
      } else {
        finalOutput = controlSubset[order(controlSubset$sampleCount, decreasing = TRUE), ][1, ]
        return(SelectionResolver(conditionDF = condition_aaa, controlDF = finalOutput))
      }
    }
    
})

# Condition Metadata
cat('Generateing Condition Metadata.\n')
matchedMetadata = as.data.frame(do.call('rbind', tempList))
head( matchedMetadata)

# Final Metadata
finalSampleMetadata = sampleMetadata
finalMatchedMetadata = matchedMetadata
  
# Variable Cleanup
#rm(case_jjj, dir_kkk)
rm(selectionVector, useDF, casesDF, controlsDF, tempList)
rm(sampleMetadata, matchedMetadata)


# Appending Condition ID
finalMatchedMetadata = data.frame(
  case_ID = paste0('Fe.', 1:nrow(finalMatchedMetadata)),
  finalMatchedMetadata
)
head(finalMatchedMetadata )
#Switching type
finalMatchedMetadata$trt_orig_sampleCount = as.numeric(finalMatchedMetadata$trt_orig_sampleCount)
finalMatchedMetadata$trt_final_sampleCount = as.numeric(finalMatchedMetadata$trt_final_sampleCount)

finalMatchedMetadata$ctl_orig_sampleCount = as.numeric(finalMatchedMetadata$ctl_orig_sampleCount)
finalMatchedMetadata$ctl_final_sampleCount = as.numeric(finalMatchedMetadata$ctl_final_sampleCount)

#Saving
cat('Now saving...\n')

#Sample Metadata
tempPath = sprintf('%sFeDrug.sample', outDir)
function.XZSaveRDS(obj = finalSampleMetadata, file = tempPath)

#Matched Metadata
tempPath = sprintf('%sFeDrug.condition', outDir)
function.XZSaveRDS(obj = finalMatchedMetadata, file = tempPath)








#2.expression matrix#####################################################################################################################
# Load Libraries
require(cmapR)
require(data.table)
require(foreach)
require(doParallel)

# Declaring Global Variables
inPath = c('Data'='level3_beta_trt_cp_n1805898x12328.gctx')
metadataDir = './Data/Compound_Data/'

options(
  stringsAsFactors = FALSE,
  warn = 1
)
# Creation of Output Directory
outDir = 'Data/Compound_Data/'
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
# To increase the read and write speed of the results, this portion is output to the memory disk
ramDir = 'ramdisk/'
dir.create(ramDir, recursive = TRUE, showWarnings = FALSE)

cat('Loading Metadata.\n')#=================================================================================
tempPath = sprintf('%sFeDrug.sample', metadataDir)
sampleDF = readRDS(tempPath)
head(sampleDF)
tempPath = sprintf('%sFeDrug.condition', metadataDir)
matchedDF = readRDS(tempPath)

# Sample Metadata Trimming
uniqueCVals = unique(c(matchedDF$trt_cval, matchedDF$ctl_cval))
sampleDF = subset(sampleDF, pert_cval %in% uniqueCVals)
rm(uniqueCVals)
#得到有用的sample
# Enable data.table optimization
sampleDF = as.data.table(sampleDF)
setkey(sampleDF, rna_centre, pert_type, pert_cval)
head(sampleDF)
dim(sampleDF)

cat('Loading Data.\n')#=======================================================================================
referenceRow = cmapR::read_gctx_ids(gctx_path = "level3_beta_trt_cp_n1805898x12328.gctx")
dataList = lapply("level3_beta_trt_cp_n1805898x12328.gctx", function(tempPath) {
  tempMatrix = parse_gctx(fname = tempPath)@mat[referenceRow, ]
  return(tempMatrix)
})
gctx_trt = dataList[[1]]
gctx_ctl<- parse_gctx("level3_beta_ctl_n188708x12328.gctx")
gctx_ctl=gctx_ctl@mat
gctxMatrix=cbind(gctx_trt,gctx_ctl)
dim(gctxMatrix)
gctxMatrix = gctxMatrix[, sampleDF$Num]



cat('Calculating...\n')# =========================================================================
# Initiate Multi-Thread
.startTime = date()
SPID_mmm = function.getPID()
referenceColumn = matchedDF$case_ID
caseCount = length(matchedDF$case_ID)
geneCount = length(referenceRow)
registerDoParallel(cores = 40)

# Parallel Processing
cat('Begin Calculating Job: ')
tempOutput = foreach(i = 1:caseCount , .inorder = TRUE) %dopar% {
  # Status Update
  if (i %% 500 == 0) {
    cat(i, ' . ', sep = '')
  }
  
  # Obtain Matching-Sample Mapping
  task_nnn = matchedDF[i, ]
  trtSamples = sampleDF[.(task_nnn$rna_centre, task_nnn$trt_type, task_nnn$trt_cval), Num]
  ctlSamples = sampleDF[.(task_nnn$rna_centre, task_nnn$ctl_type, task_nnn$ctl_cval), Num]
  
  # Sample-Size Control
  if (task_nnn$trt_orig_sampleCount != task_nnn$trt_final_sampleCount) {
    tempLogical = sample(x = 1:task_nnn$trt_orig_sampleCount, size = task_nnn$trt_orig_sampleCount, replace = FALSE)
    tempLogical = (tempLogical %in% 1:task_nnn$trt_final_sampleCount)
    trtSamples = trtSamples[tempLogical]
    rm(tempLogical)
  }
  
  if (task_nnn$ctl_orig_sampleCount != task_nnn$ctl_final_sampleCount) {
    tempLogical = sample(x = 1:task_nnn$ctl_orig_sampleCount, size = task_nnn$ctl_orig_sampleCount, replace = FALSE)
    tempLogical = (tempLogical %in% 1:task_nnn$ctl_final_sampleCount)
    ctlSamples = ctlSamples[tempLogical]
    rm(tempLogical)
  }
  #做随机抽样

  trtData = gctxMatrix[, trtSamples]
  ctlData = gctxMatrix[, ctlSamples]
  
  tempFC = rowMeans(trtData) - rowMeans(ctlData)
  
  # Save File
  tempPath = .systemInfo(eachThreadDir = ramDir, spid = SPID_mmm, type = 'FeDrug.data', idx = i)
  function.XZSaveRDS(obj = tempFC, file = tempPath)
  
  # Variable Cleanup
  rm(task_nnn, trtSamples, ctlSamples)
  rm(trtData, ctlData)
  rm(tempFC)
  function.doGC()
  return(NULL)
  #it takes time to run
}

# stop the multi-threads
registerDoSEQ()

#clean-up the variable 
#rm(matchedDF, sampleDF, gctxMatrix, allPath, tempOutput)
#function.doGC()
cat('clean-up have been done.\n')

tempMatrix = foreach(i =1:caseCount, .inorder = TRUE, .combine = cbind, .maxcombine = 1000) %do% {
  tempPath = .systemInfo(eachThreadDir = ramDir, spid = SPID_mmm, type = 'FeDrug.data', idx = i)
  return(readRDS(tempPath))
}
colnames(tempMatrix) = referenceColumn[1:caseCount]
rownames(tempMatrix) = referenceRow

tempPath = sprintf('%sFeDrug.data', outDir)
function.XZSaveRDS(obj = tempMatrix, file = tempPath)
rm(tempMatrix)
function.doGC()

# Print Timestamp
cat(sprintf('START TIME: %s\n', .startTime))
cat(sprintf('END TIME: %s\n\n', date()))

rm(matchedDF, sampleDF, gctxMatrix,  tempOutput)
################################################################################
################################################################################
FeDrug.data <- readRDS("./Data/Compound_Data/FeDrug.data")
FeDrug.data.gtc<- new("GCT", mat=FeDrug.data)
write_gctx(FeDrug.data.gtc, 
           compression_level = 9,
           "./Data/Compound_Data/FeDrug.data.gctx",
           appenddim = FALSE)

rm(list = ls())

FeDrug.condition <- readRDS("./Data/Compound_Data/FeDrug.condition")
data_trt_cp <- list()
data_trt_cp.df <- list()
fgsea.sam.trt_cp <- list()
fgsea.res.trt_cp <- list()
trt_cp_number.group <- list()
trt_cp_number <- which(FeDrug.condition$trt_type == "trt_cp")
length(trt_cp_number )
for(i in 1:56){trt_cp_number.group[[i]] <- trt_cp_number[((i-1)*10000+1):(i*10000)]}
trt_cp_number.group[[57]] <- trt_cp_number[560001:length(trt_cp_number)]
template <- parse_gctx("Data/Compound_Data/FeDrug.data.gctx",
                       rid=1:12328, cid=1:10)@mat %>% as.data.frame()
template <- template %>% dplyr::mutate(gene_id=rownames(template ))
gene_df <- fread("geneinfo_beta.txt")
template$gene_id <- as.integer(template$gene_id)
template <- dplyr::left_join(template,gene_df,by= "gene_id")
saveRDS(template,"template.RDS")
save.image("template.RData")
#========================================================================================
options(
  stringsAsFactors = FALSE,
  warn = 0
)
library(fgsea)
library(cmapR)
library(tidyverse)
resultDir="./result/cp/"
dir.create(resultDir, recursive = TRUE, showWarnings = FALSE)
load("./template.RData")
library(furrr)
library(future)
plan(multisession, workers = 40)
.startTime = date()
load("Genelist_use.Rdata")
super.third.human.pd1.all <- Genelist_use
saveRDS(super.third.human.pd1.all,"./super.third.human.pd1.all.RDS")

for(j in c(1:57)){
  setwd("workingDirectory/")
  load("./template.RData")
  fun_cmap_fgsea <- function(x){
    x <- x %>% unlist() 
    names(x) <-template$gene_symbol
    fgsea.res <- fgsea(pathways = genesets, stats = x,eps= 0.0, minSize  = 5, maxSize  = 500)
    return(fgsea.res)
  }
  genesets <- readRDS("super.third.human.pd1.all.RDS")
  data_trt_cp[[j]] <- parse_gctx("./Data/Compound_Data/FeDrug.data.gctx", 
                                 rid=1:12328, cid=trt_cp_number.group[[j]])
  data_trt_cp.df[[j]] <- as.data.frame(data_trt_cp[[j]]@mat)
  fgsea.res.trt_cp[[j]] <- furrr::future_map(data_trt_cp.df[[j]], ~ fun_cmap_fgsea(.x))
  setwd("./result/cp/");saveRDS(fgsea.res.trt_cp[[j]],paste("c",j,"_fgsea.res.trt_cp_super.RDS",sep = ""));
}

cat(sprintf('START TIME: %s\n', .startTime))
cat(sprintf('END TIME: %s\n\n', date()))
################################################################################
################################################################################
#tidy the results
setwd("workingDirectory/")
load("./template.RData")
setwd("./result/cp")

fgsea.res.trt_cp <- list()
for(i in 1:57){fgsea.res.trt_cp[[i]] <- readRDS(paste("./c",i,"_fgsea.res.trt_cp_super.RDS",sep = ""))}
fgsea.res.tidy <- list()
for(i in 1:5){print(i);
  fun_tidy= function(GSEA_res_list){
    x <- furrr::future_map_dfr(GSEA_res_list, ~ .x[i,])
    return(x)
  }
  fgsea.res.tidy[[i]] <- furrr::future_map(fgsea.res.trt_cp,~ fun_tidy(.x)) %>% Reduce(rbind,.)
}
fgsea.clue.order_with_id <- list()
for(i in 1:5){fgsea.clue.order_with_id[[i]] <- fgsea.res.tidy[[i]] %>% dplyr::mutate(id=rownames(fgsea.res.tidy[[i]]))}
for(i in 1:5){
  names(fgsea.clue.order_with_id)[i] <- fgsea.clue.order_with_id[[i]][1,1]
}
saveRDS(fgsea.clue.order_with_id,"fgsea.clue.order_with_id.RDS")
fgsea.clue.order <- readRDS("./result/cp/fgsea.clue.order_with_id.RDS")










