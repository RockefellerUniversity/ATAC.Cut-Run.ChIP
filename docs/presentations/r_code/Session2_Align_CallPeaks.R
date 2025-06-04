params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides != "yes"){
  cat("# CUT&RUN/ATAC (part 2) - Alignment and Peak Calling

---
"    
  )
  
}



## ----setwd_introtoR,eval=F----------------------------------------------------
# setwd("~/Downloads/ATAC.Cut-Run.ChIP-master/r_course")


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Aligning CUT&RUN and ATACseq reads

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Aligning CUT&RUN and ATACseq reads

---
"    
  )
  
}



## ----fa1, echo=TRUE, message = F, warning=F-----------------------------------
library(BSgenome.Mmusculus.UCSC.mm10)
BSgenome.Mmusculus.UCSC.mm10


## ----fa2,cache=FALSE,echo=TRUE------------------------------------------------
mainChromosomes <- paste0("chr",c(1:19,"X","Y","M"))
mainChrSeq <- lapply(mainChromosomes,
                     function(x)BSgenome.Mmusculus.UCSC.mm10[[x]])
names(mainChrSeq) <- mainChromosomes
mainChrSeqSet <- DNAStringSet(mainChrSeq)
mainChrSeqSet


## ----fa3, echo=TRUE,eval=FALSE------------------------------------------------
# writeXStringSet(mainChrSeqSet,
#                 "BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa")


## ----echo=TRUE,eval=FALSE-----------------------------------------------------
# library(Rsubread)
# buildindex("mm10_mainchrs","BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa",
#            memory=8000,
#            indexSplit=TRUE)
# 


## ----echo=F,eval=FALSE, include=F---------------------------------------------
# 
# chr18Seq <- BSgenome.Mmusculus.UCSC.mm10[["chr18"]]
# chr18SeqSet <- DNAStringSet(chr18Seq)
# writeXStringSet(chr18SeqSet,
#                 "BSgenome.Mmusculus.UCSC.mm10.chr18.fa")
# buildindex("mm10_chr18","BSgenome.Mmusculus.UCSC.mm10.chr18.fa",
#            memory=8000,
#            indexSplit=TRUE)


## ----eval=FALSE,include=FALSE, echo=FALSE, message = F------------------------
# require(ShortRead)
# read1 <- readFastq("~/Desktop/BRC/training/ATAC.Cut-Run.ChIP/pipeline_files/CnR/SOX9CNR_W6_rep1_allFiles/FQ_QC/SOX9CNR_W6_rep1_QC_R1.fastq.gz")
# read2 <- readFastq("~/Desktop/BRC/training/ATAC.Cut-Run.ChIP/pipeline_files/CnR/SOX9CNR_W6_rep1_allFiles/FQ_QC/SOX9CNR_W6_rep1_QC_R2.fastq.gz")
# 
# writeFastq(read1[1:1000,],"data/SOX9CNR_W6_rep1_1K_R1.fastq.gz")
# writeFastq(read2[1:1000,],"data/SOX9CNR_W6_rep1_1K_R2.fastq.gz")
# # id(read2[1:1000,])
# # myRes <- bamQC("~/Downloads/Sorted_ATAC_50K_2.bam")


## ----eval=TRUE, message = F---------------------------------------------------
require(ShortRead)

# first 1000 reads in each file
read1 <- readFastq("data/SOX9CNR_W6_rep1_1K_R1.fastq.gz")
read2 <- readFastq("data/SOX9CNR_W6_rep1_1K_R2.fastq.gz")
id(read1)[1:2]
id(read2)[1:2]


## -----------------------------------------------------------------------------
read1_toAlign <- "~/Downloads/SOX9CNR_W6_rep1_QC_R1.fastq.gz"
read2_toAlign <- "~/Downloads/SOX9CNR_W6_rep1_QC_R2.fastq.gz"


## ----echo=F,eval=TRUE, warning=F----------------------------------------------
library(Rsubread)


## ----align, echo=TRUE,eval=F, warning=F---------------------------------------
# myMapped <- align("mm10_mainchrs",
#                   readfile1 = read1_toAlign,
#                   readfile2 = read2_toAlign,
#                   type = "dna",
#                   output_file = "SOX9CNR_W6_rep1.bam",
#                   nthreads = 4,
#                   minFragLength = 0, maxFragLength = 2000)
# 


## ----sortindex, echo=TRUE,eval=FALSE, message=FALSE---------------------------
# library(Rsamtools)
# 
# sortBam("SOX9CNR_W6_rep1.bam", "SOX9CNR_W6_rep1_sorted")
# indexBam("SOX9CNR_W6_rep1_sorted.bam")


## ----coverage, echo=TRUE,eval=FALSE-------------------------------------------
# forBigWig <- coverage("SOX9CNR_W6_rep1_sorted.bam")
# forBigWig


## ----bw, echo=TRUE,eval=FALSE, message=FALSE----------------------------------
# library(rtracklayer)
# export.bw(forBigWig,con="SOX9CNR_W6_rep1.bw")


## ----weightedCover1, echo=TRUE,eval=FALSE-------------------------------------
# mappedReads <- idxstatsBam("SOX9CNR_W6_rep1_sorted.bam")
# TotalMapped <- sum(mappedReads[,"mapped"])
# 
# TotalMapped


## ----weightedCover2, echo=TRUE,eval=FALSE-------------------------------------
# 
# forBigWig <- coverage("SOX9CNR_W6_rep1_sorted.bam",
#                       weight = (10^6)/TotalMapped)
# forBigWig
# export.bw(forBigWig,con="SOX9CNR_W6_rep1_weighted.bw")


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Post alignment processing - CUT&RUN

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Post alignment processing - CUT&RUN

---
"    
  )
  
}



## ----eval=F-------------------------------------------------------------------
# BiocManager::install("Herper")
# library(Herper)
# 


## ----makeCondaEnv, echo=T, eval=F---------------------------------------------
# dir.create("miniconda")
# macs_paths <- install_CondaTools(tools= c("macs3", "samtools", "bedtools"), env="CnR_analysis", pathToMiniConda = "miniconda")


## ----eval=F, echo=F-----------------------------------------------------------
# tempdir2 <- function() {
#     tempDir <- tempdir()
#     if(dir.exists(tempDir)){
#       tempDir <- file.path(tempDir,"rr")
#     }
#     tempDir <- gsub("\\", "/", tempDir, fixed = TRUE)
#     tempDir
# }
# 
# myMiniconda <- file.path(tempdir2(), "Test")
# install_CondaTools(tools=c("macs3", "samtools", "bedtools"), env="CnR_analysis", pathToMiniConda = myMiniconda)
# 


## ----makeCondaEnv2, echo=T, eval=F--------------------------------------------
# macs_paths


## ----testCondaEnv, echo=T, eval=F---------------------------------------------
# Herper::with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
#                       system("samtools sort -h"))


## ----processData_readingInDatad, echo=TRUE,eval=TRUE,cache=FALSE--------------
library(GenomicAlignments)
flags=scanBamFlag(isProperPair = TRUE)



## ----processData_readingInDatas, echo=TRUE,eval=TRUE,cache=FALSE--------------
myParam=ScanBamParam(flag=flags,
                   what=c("qname","mapq","isize", "flag"))
myParam



## ----processData_readingInDataa, echo=TRUE,eval=T-----------------------------
sortedBAM <- "data/SOX9CNR_W6_rep1_chr18_sorted.bam"
cnrReads <- readGAlignmentPairs(sortedBAM,
                                 param=myParam)
cnrReads[1:2,]



## ----processData_readingInData2, echo=TRUE,eval=TRUE,cache=FALSE--------------
read1 <- GenomicAlignments::first(cnrReads)
read2 <- GenomicAlignments::second(cnrReads)
read2[1:2,]


## ----processData_readingInData3, echo=TRUE,eval=TRUE,cache=FALSE--------------
read1MapQ <- mcols(read1)$mapq
read2MapQ <- mcols(read2)$mapq
read1MapQ[1:5]


## ----processData_readingInData4, echo=TRUE,eval=TRUE,cache=FALSE--------------
read1MapQFreqs <- table(read1MapQ)
read2MapQFreqs <- table(read2MapQ)
read1MapQFreqs
read2MapQFreqs


## ----processData_readingInData5,fig.width=9,fig.height=4,  echo=TRUE,eval=TRUE,cache=FALSE----
library(ggplot2)
toPlot <- data.frame(MapQ=c(names(read1MapQFreqs),names(read2MapQFreqs)),
           Frequency=c(read1MapQFreqs,read2MapQFreqs),
           Read=c(rep("Read1",length(read1MapQFreqs)),rep("Read2",length(read2MapQFreqs))))
toPlot$MapQ <- factor(toPlot$MapQ,levels = unique(sort(as.numeric(toPlot$MapQ))))
ggplot(toPlot,aes(x=MapQ,y=Frequency,fill=MapQ))+
  geom_bar(stat="identity")+
  facet_grid(~Read)


## ----processData_extractingRead1, echo=T,eval=TRUE,cache=FALSE----------------
insertSizes <- abs(mcols(read1)$isize)
head(insertSizes)


## ----processData_plottingFrffagmentLengths, echo=TRUE,eval=TRUE,cache=FALSE----

fragLenSizes <- table(insertSizes)
fragLenSizes[1:5]



## ----processData_plottingFrdagmentLengths, echo=TRUE,eval=TRUE,cache=FALSE, fig.height=4, fig.width=8----
library(ggplot2)
toPlot <- data.frame(InsertSize=as.numeric(names(fragLenSizes)),
                            Count=as.numeric(fragLenSizes))
fragLenPlot <- ggplot(toPlot,aes(x=InsertSize,y=Count))+geom_line()
fragLenPlot+theme_bw()


## ----processData_plottingFragmentLengths24, echo=TRUE,eval=TRUE,cache=FALSE, fig.height=4.5, fig.width=8----
fragLenPlot+ 
  geom_vline(xintercept = c(120),colour="darkgreen")+theme_bw()



## ----processData_createOpenRegionBAM, echo=TRUE,eval=TRUE,cache=FALSE---------
cnrReads_filter <- cnrReads[insertSizes < 1000 & (!mcols(read1)$mapq == 0 | !mcols(read2)$mapq == 0)]



## ----readinBL, echo=T, eval=T-------------------------------------------------
library(rtracklayer)
blacklist <- "data/mm10-blacklist.v2.bed"
bl_regions <- rtracklayer::import(blacklist)
bl_regions


## ----removeBl, echo=T, eval=T-------------------------------------------------

fragment_spans <- granges(cnrReads_filter)
bl_remove <- overlapsAny(fragment_spans, bl_regions)
table(bl_remove)



## ----writeFilteredBAM, echo=T, eval=F-----------------------------------------
# cnrReads_filter_noBL <- cnrReads_filter[!bl_remove]
# cnrReads_unlist <- unlist(cnrReads_filter_noBL)
# names(cnrReads_unlist) <- mcols(cnrReads_unlist)$qname
# 
# filter_bam <- gsub("sorted.bam","filter.bam", basename(sortedBAM))
# rtracklayer::export(cnrReads_unlist, filter_bam,format = "bam")
# 


## ----fixmate, echo=T, eval=F--------------------------------------------------
# 
# forPeak_bam <- gsub("_filter.bam", "_forPeak.bam", filter_bam)
# Herper::with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
#                       {
#                         tempBam <- paste0(tempfile(), ".bam")
#                         system(paste("samtools sort", "-n", "-o", tempBam, filter_bam, sep = " "))
#                         system(paste("samtools fixmate", "-m", tempBam, forPeak_bam, sep = " "))
#                       })
# 


## ----bamPrepMacs, echo=T, eval=F----------------------------------------------
# 
# forMacs_bam <- gsub("_forPeak.bam", "_macs.bam", forPeak_bam)
# sortBam(forPeak_bam, gsub(".bam", "", forMacs_bam))
# indexBam(forMacs_bam)
# 


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Peak calling - CUT&RUN

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Peak calling - CUT&RUN

---
"    
  )
  
}



## ----callMacs, echo=TRUE,eval=F, warning=F------------------------------------
# peaks_name <- gsub("_macs.bam", "", basename(forMacs_bam))
# with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
#                       system2(command="macs3",args =c("callpeak",
#                       "-t", forMacs_bam,
#                       "-f", "BAMPE",
#                       "--outdir", ".",
#                       "-n", peaks_name)))
# 


## ----makeBedpe, echo=TRUE,eval=F, warning=F-----------------------------------
# # use the BAM with blacklist reads removed (from MACS section)
# bedpe <- gsub("\\.bam", "\\.bed", forPeak_bam)
# with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
#                 system(paste("bedtools bamtobed -bedpe -i", forPeak_bam, ">", bedpe, sep = " "))
# )
# 


## ----echo=T,eval=F, include = T, warning=F------------------------------------
# library(dplyr)
# library(tibble)
# 
# indexFa("BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa")
# seqlengths(Rsamtools::FaFile("BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa")) %>%
#   as.data.frame() %>%
#   rownames_to_column() %>%
#   write.table(file = "chrom.lengths.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


## ----makeBedgraph, echo=T,eval=F, warning=F-----------------------------------
# 
# bedgraph <- gsub("\\.bed", "\\.bedgraph", bedpe)
# with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
#                 system(paste("bedtools genomecov -bg -i", bedpe, "-g", "data/chrom.lengths.txt", ">", bedgraph, sep = " "))
# )


## ----getSEACR, echo=T,eval=F, warning=F---------------------------------------
# download.file("https://github.com/FredHutch/SEACR/archive/refs/tags/v1.4-beta.2.zip", destfile = "~/Downloads/SEACR_v1.4.zip")
# unzip("~/Downloads/SEACR_v1.4.zip", exdir = "~/Downloads/SEACR_v1.4" )
# seacr_path <- "~/Downloads/SEACR_v1.4/SEACR-1.4-beta.2/SEACR_1.4.sh"
# system(paste(seacr_path, "-h"))


## ----SEACRpermissions, echo=T,eval=F, warning=F-------------------------------
# system(paste("chmod 777", seacr_path))
# system(paste(seacr_path, "-h"))


## ----runSEACR, echo=T,eval=F, warning=F---------------------------------------
# seacr_prefix <- gsub("\\.bedgraph", "_top01", bedgraph)
# 
# system(paste(seacr_path, "-b", bedgraph, "-c 0.01", "-n non", "-m  stringent", "-o", seacr_prefix))


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Post alignment processing and peak calling - ATACseq

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Post alignment processing and peak calling - ATACseq

---
"    
  )
  
}



## ----echo=TRUE,eval=TRUE,cache=FALSE------------------------------------------

flags=scanBamFlag(isProperPair = TRUE)
myParam=ScanBamParam(flag=flags,
                   what=c("qname","mapq","isize", "flag"))

sortedBAM_atac <- "data/W6_ATAC_rep1_chr18_sorted.bam"
atacReads <- readGAlignmentPairs(sortedBAM_atac,
                                 param=myParam)
head(atacReads, 2)


## ----processData_plottingFrffagmentLengths_atac, echo=TRUE,eval=TRUE,cache=FALSE,fig.width=6,fig.height=4----
atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes_atac <- abs(elementMetadata(atacReads_read1)$isize)
fragLenSizes_atac <- table(insertSizes_atac)

toPlot_atac <- data.frame(InsertSize=as.numeric(names(fragLenSizes_atac)),
                            Count=as.numeric(fragLenSizes_atac))
fragLenPlot_atac <- ggplot(toPlot_atac,aes(x=InsertSize,y=Count))+geom_line() +theme_bw()


## ----processData_plottingFragmentLengths3_atac, echo=TRUE,eval=TRUE,cache=FALSE,fig.width=6,fig.height=4----

fragLenPlot_atac + scale_y_continuous(trans='log2')


## ----processData_plottingFragmentLengths24_atac, echo=TRUE,eval=TRUE,cache=FALSE,fig.width=6,fig.height=4----
fragLenPlot_atac+ scale_y_continuous(trans='log2')+
  geom_vline(xintercept = c(180,247),colour="red")+
  geom_vline(xintercept = c(315,437),colour="darkblue")+
  geom_vline(xintercept = c(100),colour="darkgreen")+theme_bw()



## ----processData_createOpenRegionBAM_atac, echo=TRUE,eval=TRUE,cache=FALSE----
atacReads_filter <- atacReads[!mcols(read1)$mapq == 0 | !mcols(read2)$mapq == 0]



## ----removeBl_atac, echo=T, eval=F--------------------------------------------
# # the 'bl_regions' GRanges object was generated in CUT&RUN section
# fragment_spans_atac <- granges(atacReads_filter)
# bl_remove_atac <- overlapsAny(fragment_spans_atac, bl_regions)
# table(bl_remove_atac)
# 


## ----writeFilteredBAM_atac, echo=T, eval=F------------------------------------
# atacReads_filter_noBL <- atacReads_filter[!bl_remove]
# atacReads_unlist <- unlist(atacReads_filter_noBL)
# names(atacReads_unlist) <- mcols(atacReads_unlist)$qname
# 
# filter_bam_atac <- gsub("sorted.bam","filter.bam", basename(sortedBAM_atac))
# rtracklayer::export(atacReads_unlist, filter_bam_atac,format = "bam")
# 


## ----fixmate_atac, echo=T, eval=F---------------------------------------------
# 
# forPeak_bam_atac <- gsub("_filter.bam", "_macs.bam", filter_bam_atac)
# Herper::with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
#                       {
#                         tempBam <- paste0(tempfile(), ".bam")
#                         tempBam2 <- paste0(tempfile(), ".bam")
#                         system(paste("samtools sort", "-n", "-o", tempBam, filter_bam_atac, sep = " "))
#                         system(paste("samtools fixmate", "-m", tempBam, tempBam2, sep = " "))
#                         system(paste("samtools sort", "-n", "-o", forPeak_bam_atac, tempBam2, sep = " "))
#                       })
# 


## ----callMacs_atac, echo=TRUE,eval=F, warning=F-------------------------------
# peaks_name_atac <- gsub("_macs.bam", "", basename(forPeak_bam_atac))
# with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
#                       system2(command="macs3",args =c("callpeak",
#                       "-t", forPeak_bam_atac,
#                       "--format", "BAMPE",
#                       "--outdir", ".",
#                       "--name", peaks_name_atac)))
# 


## ----writeNFBAM_atac, echo=T, eval=F------------------------------------------
# read1_noBL <- GenomicAlignments::first(atacReads_filter_noBL)
# insertSizes_atac_noBL <- abs(elementMetadata(read1_noBL)$isize)
# 
# atacReads_NF <- atacReads_filter_noBL[insertSizes_atac_noBL < 100]
# atacReads_NF_unlist <- unlist(atacReads_NF)
# names(atacReads_NF_unlist) <- mcols(atacReads_NF_unlist)$qname
# 
# NF_bam_atac <- gsub("sorted.bam","filterNF.bam", basename(sortedBAM_atac))
# rtracklayer::export(atacReads_NF_unlist, NF_bam_atac,format = "bam")
# 


## ----fixmate_atacNF, echo=T, eval=F-------------------------------------------
# 
# forPeak_NFbam_atac <- gsub("_filterNF.bam", "_macsNF.bam", NF_bam_atac)
# Herper::with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
#                       {
#                         tempBam <- paste0(tempfile(), ".bam")
#                         tempBam2 <- paste0(tempfile(), ".bam")
#                         system(paste("samtools sort", "-n", "-o", tempBam, NF_bam_atac, sep = " "))
#                         system(paste("samtools fixmate", "-m", tempBam, tempBam2, sep = " "))
#                         system(paste("samtools sort", "-n", "-o", forPeak_NFbam_atac, tempBam2, sep = " "))
#                       })
# 


## ----callMacs_atacNF, echo=TRUE,eval=F, warning=F-----------------------------
# peaksNF_name_atac <- gsub("_macsNF.bam", "_macsNF", basename(forPeak_NFbam_atac))
# with_CondaEnv("CnR_analysis", pathToMiniConda = "miniconda",
#                       system2(command="macs3",args =c("callpeak",
#                       "-t", forPeak_NFbam_atac,
#                       "--format", "BAMPE",
#                       "--outdir", ".",
#                       "--name", peaksNF_name_atac)))
# 

