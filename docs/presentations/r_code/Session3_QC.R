params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T) # delete cache before any merging 



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides != "yes"){
  cat("# CUT&RUN/ATAC (part 3) - Quality Control

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Quality control - CUT&RUN

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Quality control - CUT&RUN

---
"    
  )
  
}


## ----eval=F, echo=F-----------------------------------------------------------
# mappedReads <- idxstatsBam("../../../../pipeline_files/CnR/BAMS_forChIPQC/SOX9CNR_W6_rep1_sorted.bam")
# mappedReads <- mappedReads[mappedReads$seqnames %in% paste0("chr", c(1:19, "X", "Y", "M")), ]
# saveRDS(mappedReads, file="data/idxstatsBam_sox9_CnR_W6R1.rds")


## ----mapped1, echo=TRUE,eval=FALSE--------------------------------------------
# mappedReads <- idxstatsBam("~/Downloads/SOX9CNR_W6_rep1_sorted.bam")
# mappedReads <- mappedReads[mappedReads$seqnames %in% paste0("chr", c(1:19, "X", "Y", "M")), ]
# TotalMapped <- sum(mappedReads[,"mapped"])
# ggplot(mappedReads,aes(x=seqnames,y=mapped))+
#   geom_bar(stat="identity")+coord_flip()


## ----mapped, echo=FALSE,eval=TRUE,fig.width=4,fig.height=4--------------------
mappedReads <- readRDS("data/idxstatsBam_sox9_CnR_W6R1.rds")
TotalMapped <- sum(mappedReads[,"mapped"])
suppressPackageStartupMessages(library(ggplot2))
ggplot(mappedReads,aes(x=seqnames,y=mapped))+geom_bar(stat="identity")+coord_flip()


## ----mycQCdwdwshowL,include=FALSE---------------------------------------------
library(ChIPQC)


## ----eval=F-------------------------------------------------------------------
# QCresult <- ChIPQCsample(reads="/pathTo/myCnRreads.bam",
#                          genome="mm10",
#                          peaks = "/pathTo/myCnRpeaks.bed",
#                          blacklist = "/pathTo/mm10_Blacklist.bed")


## ----mycQC,cache=F,eval=TRUE, message=F---------------------------------------
library(ChIPQC)
blklist <- rtracklayer::import.bed("data/mm10-blacklist.v2.bed")
qc_sox9_rep1 <- ChIPQCsample("data/SOX9CNR_W6_rep1_chr18_sorted.bam",
                             annotation = "mm10",
                             peaks = "data/SOX9CNR_W6_rep1_chr18_peaks.narrowPeak",
                             blacklist = blklist,
                             chromosomes = "chr18")


## ----mycQC2,cache=F,eval=TRUE-------------------------------------------------
class(qc_sox9_rep1)


## ----mycQCshow,eval=TRUE------------------------------------------------------
qc_sox9_rep1


## ----mycQCshow2,cache=F,eval=FALSE, echo=T------------------------------------
# bamsToQC <- c("~/Downloads/SOX9CNR_D0_rep1_sorted.bam",
#               "~/Downloads/SOX9CNR_D0_rep2_sorted.bam",
#               "~/Downloads/SOX9CNR_W6_rep1_sorted.bam",
#               "~/Downloads/SOX9CNR_W6_rep2_sorted.bam")
# 
# peaksToQC <- c("data/peaks/SOX9CNR_D0_rep1_macs_peaks.narrowPeak",
#               "data/peaks/SOX9CNR_D0_rep2_macs_peaks.narrowPeak",
#               "data/peaks/SOX9CNR_W6_rep1_macs_peaks.narrowPeak",
#               "data/peaks/SOX9CNR_W6_rep2_macs_peaks.narrowPeak")
# 


## ----mycQCshow4,cache=F,eval=FALSE, echo=T------------------------------------
# 
# myQC <- lapply(seq_along(bamsToQC),function(x){
#   ChIPQCsample(
#     bamsToQC[x],
#     annotation = "mm10",
#     peaks = peaksToQC[x],
#     blacklist = blklist,
#     chromosomes = "chr18"
#   )
# })
# names(myQC) <- basename(bamsToQC)
# 


## ----mycQCshow3,cache=F,eval=FALSE, echo=F------------------------------------
# bamsToQC <- c("../../../../pipeline_files/CnR/BAMS_forChIPQC/SOX9CNR_D0_rep1_sorted.bam",
#               "../../../../pipeline_files/CnR/BAMS_forChIPQC/SOX9CNR_D0_rep2_sorted.bam",
#               "../../../../pipeline_files/CnR/BAMS_forChIPQC/SOX9CNR_W6_rep1_sorted.bam",
#               "../../../../pipeline_files/CnR/BAMS_forChIPQC/SOX9CNR_W6_rep2_sorted.bam")
# 
# peaksToQC <- c("data/peaks/SOX9CNR_D0_rep1_macs_peaks.narrowPeak",
#               "data/peaks/SOX9CNR_D0_rep2_macs_peaks.narrowPeak",
#               "data/peaks/SOX9CNR_W6_rep1_macs_peaks.narrowPeak",
#               "data/peaks/SOX9CNR_W6_rep2_macs_peaks.narrowPeak")
# 
# 
# myQC <- lapply(seq_along(bamsToQC),function(x){
#   ChIPQCsample(
#     bamsToQC[x],
#     annotation = "mm10",
#     peaks = peaksToQC[x],
#     blacklist = blklist,
#     chromosomes = "chr18"
#   )
# })
# names(myQC) <- basename(bamsToQC)
# 
# saveRDS(myQC, "sox9_QC_withPeaks.rds")


## ----qcmetricsA,include=FALSE, echo=F, eval=T---------------------------------
myQC <- readRDS("data/sox9_QC_withPeaks.rds")


## ----qcmetrics,cache=FALSE,eval=TRUE------------------------------------------
QCmetrics(myQC)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Blacklists and SSD

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Blacklists and SSD

---
"    
  )
  
}



## ----fig.width=6,fig.height=2,warning=FALSE,message=FALSE---------------------
plotSSD(myQC)+xlim(0,7)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Library complexity and enrichment

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Library complexity and enrichment

---
"    
  )
  
}



## ----fig.width=6,fig.height=3,warning=FALSE,message=FALSE---------------------
myFlags <- flagtagcounts(myQC)
myFlags["DuplicateByChIPQC",]/myFlags["Mapped",]


## ----warning=FALSE,message=FALSE,fig.width=8,fig.height=4---------------------
p <- plotRegi(myQC)


## ----warning=FALSE,fig.width=12,fig.height=6----------------------------------
p


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Quality control - ATACseq

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Quality control - ATACseq

---
"    
  )
  
}


## ----eval=F, echo=F-----------------------------------------------------------
# mappedReads_atac <- idxstatsBam("../../../../pipeline_files/ATAC/W6_ATAC_rep1_allFiles/BAM/W6_ATAC_rep1_sorted.bam")
# mappedReads_atac <- mappedReads_atac[mappedReads_atac$seqnames %in% paste0("chr", c(1:19, "X", "Y", "M")), ]
# saveRDS(mappedReads_atac,file="data/ATAC_WTR1_idxstats.rds")


## ----eval=T, echo=F,cache=F,dependson="processData_setBAM"--------------------
mappedReads_atac <- readRDS("data/ATAC_WTR1_idxstats.rds")


## ----quickMappingStatsPerChromosomea, echo=TRUE,eval=F------------------------
# library(Rsamtools)
# mappedReads_atac <- idxstatsBam("W6_ATAC_rep1_sorted.bam")
# mappedReads_atac <- mappedReads_atac[mappedReads_atac$seqnames %in% paste0("chr", c(1:19, "X", "Y", "M")), ]
# 


## ----quickMappingStatsPerChromosomes, echo=TRUE,eval=TRUE,cache=F,dependson="processData_setBAM",fig.height=4,fig.width=6----
library(ggplot2)

ggplot(mappedReads_atac,aes(seqnames,mapped))+
  geom_bar(stat="identity")+coord_flip()


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# ATACseqQC

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# ATACseqQC

---
"    
  )
  
}


## ----eval=FALSE---------------------------------------------------------------
# BiocManager::install("ATACseqQC")


## ----echo=T, eval=T-----------------------------------------------------------
library(ATACseqQC)

atac_bam <- "data/W6_ATAC_rep1_chr18_sorted.bam"
ATACQC <- bamQC(atac_bam)


## -----------------------------------------------------------------------------
names(ATACQC)


## -----------------------------------------------------------------------------
ATACQC$PCRbottleneckCoefficient_1


## -----------------------------------------------------------------------------
ATACQC$PCRbottleneckCoefficient_2


## ----fig.height=6, fig.width=6------------------------------------------------

fragSize <- fragSizeDist(atac_bam,
                         bamFiles.labels = gsub("\\.bam", "", basename(atac_bam)))


## -----------------------------------------------------------------------------
fragSize$W6_ATAC_rep1_chr18_sorted[1:10]


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Evaluating TSS signal

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Evaluating TSS signal

---
"    
  )
  
}


## ----processData_aligna, echo=TRUE,eval=FALSE,cache=FALSE---------------------
# BiocManager::install("soGGi")
# library(soGGi)


## ----processData_txdb, echo=TRUE,eval=TRUE,cache=FALSE, message=F-------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
TxDb.Mmusculus.UCSC.mm10.knownGene


## ----processData_genes, echo=TRUE,eval=TRUE,cache=FALSE-----------------------
genesLocations <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
genesLocations


## ----processData_resize, echo=TRUE,eval=TRUE,cache=FALSE----------------------
tssLocations <- resize(genesLocations,fix="start",width = 1)
tssLocations



## ----processData_subset, echo=TRUE,eval=T-------------------------------------
mainChromosomes <- paste0("chr",c(1:19,"X","Y","M"))

tssLocations <- tssLocations[as.vector(seqnames(tssLocations)) %in% mainChromosomes]

seqlevels(tssLocations) <- mainChromosomes



## ----processData_soggi, echo=TRUE,eval=FALSE,cache=FALSE----------------------
# library(soGGi)
# sortedBAM <- "~/Downloads/W6_ATAC_rep1_sorted.bam"
# 
# library(Rsamtools)
# # Nucleosome free
# allSignal <- regionPlot(bamFile = sortedBAM,
#                         testRanges = tssLocations)


## ----processData_soggia, echo=TRUE,eval=FALSE,cache=FALSE---------------------
# nucFree <- regionPlot(bamFile = sortedBAM,
#                         testRanges = tssLocations,
#                         style = "point",
#                         format="bam",
#                         paired=TRUE,
#                         minFragmentLength = 0,
#                         maxFragmentLength = 100,
#                         forceFragment = 50)
# 


## ----processData_soggia1, echo=F,eval=F---------------------------------------
# saveRDS(nucFree, "data/nucFree_TSS.rds")


## ----processData_soggia2, echo=F,eval=T---------------------------------------
library(soGGi)
nucFree <- readRDS("data/nucFree_TSS.rds")


## ----processData_plot,fig.height=3.5,fig.width=7, echo=TRUE,eval=T,cache=F,message=FALSE,warning=FALSE----
plotRegion(nucFree)


## ----processData_soggi4, echo=F,eval=FALSE,cache=FALSE,message=FALSE,warning=FALSE----
# monoNuc <- regionPlot(bamFile = sortedBAM,
#                         testRanges = tssLocations,
#                         style = "point",
#                         format="bam",
#                         paired=TRUE,
#                         minFragmentLength = 180,maxFragmentLength = 240,forceFragment = 80)
# saveRDS(monoNuc, "data/monoNuc_TSS.rds")


## ----processData_soggi2.5, echo=F,eval=T,cache=FALSE,message=FALSE,warning=FALSE----

monoNuc <- readRDS(file = "data/monoNuc_TSS.rds")


## ----processData_soggi3, echo=TRUE,eval=FALSE,cache=FALSE,message=FALSE,warning=FALSE----
# monoNuc <- regionPlot(bamFile = sortedBAM,
#                         testRanges = tssLocations,
#                         style = "point",
#                         format="bam",
#                         paired=TRUE,
#                         minFragmentLength = 180,maxFragmentLength = 240,forceFragment = 80)
# 


## ----processData_plot3, fig.height=3.5,fig.width=7, echo=TRUE,eval=T,cache=FALSE, message=FALSE, warning=FALSE----
plotRegion(monoNuc)

