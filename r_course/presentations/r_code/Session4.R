params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides != "yes"){
  cat("# Epigenomics (part 4)

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Consensus Peaks

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Consensus Peaks

---
"    
  )
  
}



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
peakFiles <- dir("data/peaks",pattern="narrowPeak",
                 full.names = TRUE)
peakFiles



## ----eval=T,echo=T, warning=FALSE---------------------------------------------
library(rtracklayer)
macsPeaks_list <- lapply(peakFiles, import)
length(macsPeaks_list)


## -----------------------------------------------------------------------------
macsPeaks_list[[1]]



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
macsPeaks_GRL <- GRangesList(macsPeaks_list)
names(macsPeaks_GRL) <- c("D0_rep1","D0_rep2","W6_rep1","W6_rep2")
class(macsPeaks_GRL)
names(macsPeaks_GRL)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
lengths(macsPeaks_GRL)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
macsPeaks_GRLCentred <- resize(macsPeaks_GRL,10,fix="center")
width(macsPeaks_GRLCentred)


## ----eval=T,echo=T, warning=FALSE---------------------------------------------
D0_rep1_Peaks <- macsPeaks_GRL$D0_rep1
D0_rep2_Peaks <- macsPeaks_GRL$D0_rep2
length(D0_rep1_Peaks)
length(D0_rep2_Peaks)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
D0_rep1_unique.Peaks <- D0_rep1_Peaks[!D0_rep1_Peaks %over% D0_rep2_Peaks]
D0_rep2_unique.Peaks <- D0_rep2_Peaks[!D0_rep2_Peaks %over% D0_rep1_Peaks]
length(D0_rep1_unique.Peaks)
length(D0_rep2_unique.Peaks)
export.bed(D0_rep1_unique.Peaks,"D0_rep1_Unique.bed")
export.bed(D0_rep2_unique.Peaks,"D0_rep2_Unique.bed")


## ----eval=T,echo=T, warning=FALSE---------------------------------------------

D0_rep1_common.Peaks <- D0_rep1_Peaks[D0_rep1_Peaks %over% D0_rep2_Peaks]
D0_rep2_common.Peaks <- D0_rep2_Peaks[D0_rep2_Peaks %over% D0_rep1_Peaks]
length(D0_rep1_common.Peaks)
length(D0_rep2_common.Peaks)
export.bed(D0_rep1_common.Peaks,"D0_rep1_Common.bed")
export.bed(D0_rep2_common.Peaks,"D0_rep2_Common.bed")




## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
allPeaksSet_Overlapping <- unlist(macsPeaks_GRL)
allPeaksSet_Overlapping


## ----eval=T,echo=T, warning=FALSE---------------------------------------------
allPeaksSet_nR <- reduce(allPeaksSet_Overlapping)
allPeaksSet_nR
export.bed(allPeaksSet_nR,"allPeaksSet_nR.bed")


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
commonPeaks <- allPeaksSet_nR[allPeaksSet_nR %over% D0_rep1_Peaks &
                               allPeaksSet_nR %over% D0_rep2_Peaks]
commonPeaks
export.bed(commonPeaks,"D0_commonPeaks.bed")



## ----eval=T,echo=T,warning=FALSE----------------------------------------------
overlap <- list()
for(i in 1:length(macsPeaks_GRL)){
  overlap[[i]] <- allPeaksSet_nR %over% macsPeaks_GRL[[i]]
}
overlap[[1]][1:2]


## ----eval=T,echo=T, warning=FALSE---------------------------------------------
overlapMatrix <- do.call(cbind,overlap)
colnames(overlapMatrix) <- names(macsPeaks_GRL)
overlapMatrix[1:3,]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
mcols(allPeaksSet_nR) <- overlapMatrix
allPeaksSet_nR[1:3,]


## ----eval=T,echo=T, warning=FALSE,fig.height=5,fig.width=5--------------------
library(limma)
vennDiagram(mcols(allPeaksSet_nR))


## ----eval=T,echo=T, warning=FALSE---------------------------------------------
vennCounts(mcols(allPeaksSet_nR))


## ----eval=T,echo=T,warning=FALSE----------------------------------------------
D0_HC_Peaks <- allPeaksSet_nR[rowSums(as.data.frame(mcols(allPeaksSet_nR)[,c("D0_rep1","D0_rep2")])) >= 2]

tail(D0_HC_Peaks)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
D0_HC_UniquePeaks <- allPeaksSet_nR[
  rowSums(as.data.frame(
    mcols(allPeaksSet_nR)[,c("D0_rep1","D0_rep2")])) >= 2 &
  rowSums(as.data.frame(
    mcols(allPeaksSet_nR)[,c("W6_rep1","W6_rep2")])) == 0  
  ]
export.bed(D0_HC_UniquePeaks,"D0_HC_UniquePeaks.bed")
D0_HC_UniquePeaks[1:5,]


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Differential Peaks

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Differential Peaks

---
"    
  )
  
}



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
HC_Peaks <- allPeaksSet_nR[
  rowSums(as.data.frame(
    mcols(allPeaksSet_nR)[,c("D0_rep1","D0_rep2")])) >= 2 |
  rowSums(as.data.frame(
    mcols(allPeaksSet_nR)[,c("W6_rep1","W6_rep2")])) >= 2  
  ]
export.bed(HC_Peaks,"HC_Peaks.bed")
HC_Peaks


## ----eval=F, echo=F, warning=FALSE--------------------------------------------
# 
# library(Rsamtools)
# 
# datasets <- rbind(
# c("SOX9CNR_D0_rep1",	"YHBLDAMZIWFJXNRSVKOT"),
# c("SOX9CNR_D0_rep2",	"ZYCLMWGQSVOHTNAIFPDX"),
#   c("SOX9CNR_W6_rep1",	"LIQKNMDFWPVBJSTHOEYG"),
#     c("SOX9CNR_W6_rep2",	"TWHECASNYQDIKUOXPVGB"))
# 
# bams <- sapply(1:4, function(x){
# 
# path <- dir(file.path("/rugpfs/fs0/brc/scratch/brc_pipeline/analysis/testBRC_Shiny/",datasets[x,2],"workflow_data","BAM"), pattern=".bam$", full.names=T)
# 
# return(path)
# })
# 
# bamFL <- BamFileList(bams,yieldSize = 5000000)
# bamFL


## ----eval=F, echo=T, warning=FALSE--------------------------------------------
# 
# library(Rsamtools)
# 
# bams <- c("SOX9CNR_D0_rep1.fastq.gz",
#           "SOX9CNR_D0_rep2.fastq.gz",
#           "SOX9CNR_W6_rep1.fastq.gz",
#           "SOX9CNR_W6_rep2.fastq.gz")
# 
# bamFL <- BamFileList(bams,yieldSize = 5000000)
# bamFL


## ----eval=F, echo=F, warning=FALSE--------------------------------------------
# HC_Peaks <- rtracklayer::import("HC_Peaks.bed")
# library(GenomicAlignments)
# MyCounts <- summarizeOverlaps(HC_Peaks,
#                               reads = bamFL,
#                               ignore.strand = TRUE)
# save(MyCounts,file="MyCounts.RData")


## ----eval=F, echo=T, warning=FALSE--------------------------------------------
# library(GenomicAlignments)
# MyCounts <- summarizeOverlaps(HC_Peaks,
#                               reads = bamFL,
#                               ignore.strand = TRUE)
# save(MyCounts,file="MyCounts.RData")


## ----warning=FALSE------------------------------------------------------------
load("data/MyCounts.RData")

class(MyCounts)

MyCounts



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
metaDataFrame <- data.frame(DevStage=c("D0","D0","W6","W6"))
rownames(metaDataFrame) <- colnames(MyCounts)
metaDataFrame


## ----eval=T,echo=T, warning=FALSE---------------------------------------------
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = assay(MyCounts),
                              colData = metaDataFrame,
                              design = ~ DevStage,
                              rowRanges= HC_Peaks)


## ----eval=T,echo=T, warning=FALSE---------------------------------------------
dds <- DESeq(dds)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
dds


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
W6MinusD0 <- results(dds,
                        contrast = c("DevStage","W6","D0"),
                        format="GRanges")
W6MinusD0 <- W6MinusD0[order(W6MinusD0$pvalue),]
W6MinusD0


## ----eval=T,echo=T, warning=FALSE---------------------------------------------
W6MinusD0


## -----------------------------------------------------------------------------
rio::export(as.data.frame(W6MinusD0), "data/W6MinusD0.xlsx")




## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
W6MinusD0.Filt <- W6MinusD0[!is.na(W6MinusD0$pvalue) | !is.na(W6MinusD0$padj)]
UpinW6 <- W6MinusD0[W6MinusD0$padj < 0.05 & W6MinusD0$log2FoldChange > 0]
DowninW6 <- W6MinusD0[W6MinusD0$padj < 0.05 & W6MinusD0$log2FoldChange < 0]
export.bed(UpinW6,"UpinW6.bed")
export.bed(DowninW6,"DowninW6.bed")



## ----eval=T,echo=T, eval=F, echo=T, warning=FALSE-----------------------------
# library(tracktables)
# myReport <- makebedtable(W6MinusD0.Filt,"W6MinusD0.html",
#                          basedirectory = getwd())
# 
# browseURL(myReport)


## -----------------------------------------------------------------------------
myrlog <- rlog(dds)
plotPCA(myrlog, intgroup="DevStage")


## -----------------------------------------------------------------------------
library(EnhancedVolcano)

EnhancedVolcano(as.data.frame(W6MinusD0.Filt),
                lab=paste0(seqnames(W6MinusD0.Filt), ranges(W6MinusD0.Filt)),
                x = "log2FoldChange",
                y="padj", selectLab = "")



