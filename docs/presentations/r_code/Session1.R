params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides != "yes"){
  cat("# Epigenomics

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Set Up

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Set Up

---
"    
  )
  
}



## ----setwd_introtoR,eval=F----------------------------------------------------
# setwd("/PathToMyDownload/RU_Course_template/r_course")
# # e.g. setwd("~/Downloads/Intro_To_R_1Day/r_course")


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Epigenomic Approaches

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Epigenomic Approaches

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Fastq QC

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Fastq QC

---
"    
  )
  
}



## ----shortreada,include=FALSE-------------------------------------------------
library(ShortRead)



## ----shortread, warning=F, message=F------------------------------------------
library(ShortRead)


## ----echo=F,eval=F------------------------------------------------------------
# fqSample <- FastqSampler("~/Downloads/SOX9CNR_D0_rep1_R1.fastq.gz",n=10^6)
# fastq <- yield(fqSample)
# 
# writeFastq(fastq,file = "../data/SOX9CNR_D0_rep1_R1_subsample.fastq.gz",mode = "w")
# 


## ----eval=F, echo=F-----------------------------------------------------------
# fastq <- readFastq(dirPath = "data/SOX9CNR_D0_rep1_R1_subsample.fastq.gz")


## ----mycRep1Reads,echo=T,eval=F-----------------------------------------------
# fqSample <- FastqSampler("~/Downloads/SOX9CNR_D0_rep1_R1.fastq.gz",n=10^6)
# fastq <- yield(fqSample)


## ----eval=T, echo=T-----------------------------------------------------------
fastq <- readFastq(dirPath = "data/SOX9CNR_D0_rep1_R1_subsample.fastq.gz")


## ----mycRep1ReadsShortReadQ---------------------------------------------------
fastq


## ----mycRep1ReadsAccessor-----------------------------------------------------
readSequences <- sread(fastq)
readQuality <- quality(fastq)
readIDs <- id(fastq)
readSequences


## ----mycRep1ReadsQScores------------------------------------------------------
readQuality <- quality(fastq)
readQualities <- alphabetScore(readQuality)
readQualities[1:10]


## ----mycRep1ReadsQScoresPlot--------------------------------------------------
library(ggplot2)
toPlot <- data.frame(ReadQ=readQualities)
ggplot(toPlot,aes(x=ReadQ))+geom_histogram()+theme_minimal()


## ----mycRep1ReadsAlpFreq------------------------------------------------------
readSequences <- sread(fastq)
readSequences_AlpFreq <- alphabetFrequency(readSequences)
readSequences_AlpFreq[1:3,]


## ----mycRep1ReadsAlpFreqSum---------------------------------------------------
summed__AlpFreq  <- colSums(readSequences_AlpFreq)
summed__AlpFreq[c("A","C","G","T","N")]


## ----mycRep1ReadsAlpByCycle---------------------------------------------------
readSequences_AlpbyCycle <- alphabetByCycle(readSequences)
readSequences_AlpbyCycle[1:4,1:10]


## ----mycRep1ReadsAlpByCyclePlot-----------------------------------------------
AFreq <- readSequences_AlpbyCycle["A",]
CFreq <- readSequences_AlpbyCycle["C",]
GFreq <- readSequences_AlpbyCycle["G",]
TFreq <- readSequences_AlpbyCycle["T",]
toPlot <- data.frame(Count=c(AFreq,CFreq,GFreq,TFreq),
                     Cycle=rep(1:40,4),
                     Base=rep(c("A","C","G","T"),each=40))



## ----mycRep1ReadsAlpByCyclePlot2,cache=TRUE,eval=FALSE,dependson="mycRep1ReadsAlpByCyclePlot",fig.height=4,fig.width=8----
# 
# ggplot(toPlot, aes(y=Count,x=Cycle,colour=Base)) + geom_line() +
#   theme_bw()


## ----mycRep1ReadsAlpByCyclePlot3----------------------------------------------

ggplot(toPlot,aes(y=Count,x=Cycle,colour=Base)) + geom_line() + ylim(150000,400000) +
  theme_bw()


## ----mycRep1ReadsQByCycle,cache=TRUE,dependson="mycRep1ReadsAlpFreq"----------
qualAsMatrix <- as(readQuality,"matrix")
qualAsMatrix[1:2,]


## ----mycRep1ReadsQByCyclePlot,cache=TRUE,dependson="mycRep1ReadsQByCycle",fig.width=8,fig.height=4----
boxplot(qualAsMatrix[1:1000,])


## ----eval=F-------------------------------------------------------------------
# library(Rfastp)
# 
# rfastp(read1="/Users/mattpaul/Downloads/SRR20110418_1.fastq.gz",
#        read2 = "/Users/mattpaul/Downloads/SRR20110418_2.fastq.gz",
#        outputFastq = "SOX9CNR_D0_rep1_filtered")

