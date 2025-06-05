params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
suppressPackageStartupMessages(require(knitr))
suppressPackageStartupMessages(require(memes))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides != "yes"){
  cat("# Epigenomics (part 6)

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Motif Databases

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Motif Databases

---
"    
  )
  
}


## ----eval=TRUE----------------------------------------------------------------
library(MotifDb)
MotifDb


## ----eval=TRUE----------------------------------------------------------------
class(MotifDb)


## ----eval=TRUE----------------------------------------------------------------
length(MotifDb)
MotifNames <- names(MotifDb)
MotifNames[1:10]


## ----eval=TRUE----------------------------------------------------------------
MotifDb[1]


## ----eval=TRUE----------------------------------------------------------------
MotifDb[[1]]
colSums(MotifDb[[1]])


## ----eval=TRUE----------------------------------------------------------------
values(MotifDb)[1:2,]


## ----eval=TRUE----------------------------------------------------------------
Sox9Motifs <- query(MotifDb,"Sox9")
Sox9Motifs


## ----eval=TRUE----------------------------------------------------------------
Sox9Motifs <- query(MotifDb,c("Sox9","hsapiens","jaspar2022"))
Sox9Motifs


## ----eval=TRUE----------------------------------------------------------------
library(JASPAR2024)
library(RSQLite)

JASPAR2024 <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))




## ----eval=TRUE, message=F, warning=F------------------------------------------
library(TFBSTools)
?getMatrixByID


## ----eval=TRUE,echo=TRUE------------------------------------------------------
SOX9mat <- getMatrixByName(sq24,"SOX9")
class(SOX9mat)



## ----eval=TRUE,echo=TRUE------------------------------------------------------
SOX9mat2 <- getMatrixByID(sq24,"MA0077.1")


## ----eval=TRUE,echo=TRUE------------------------------------------------------
ID(SOX9mat)


## ----eval=TRUE----------------------------------------------------------------
myMatrix <- Matrix(SOX9mat)
myMatrixToo <- as.matrix(myMatrix)
myMatrix


## ----eval=TRUE----------------------------------------------------------------
opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "vertebrates"
opts[["all_versions"]] <- FALSE
motifList <- getMatrixSet(sq24, opts)
   


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Visualizing motifs

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Visualizing motifs

---
"    
  )
  
}


## ----eval=TRUE,echo=FALSE,warning=FALSE,message=FALSE, fig.height=4,fig.width=6----
require(seqLogo)
Sox9Motifs <- query(MotifDb,"Sox9")
seqLogo::seqLogo(Sox9Motifs[[1]],ic.scale = FALSE)


## ----eval=TRUE, fig.height=4,fig.width=6--------------------------------------
library(seqLogo)
Sox9Motifs <- query(MotifDb,"Sox9")
seqLogo::seqLogo(Sox9Motifs[[1]], ic.scale = FALSE)


## ----eval=TRUE, fig.height=4,fig.width=6--------------------------------------
seqLogo::seqLogo(Sox9Motifs[[1]])


## ----eval=TRUE----------------------------------------------------------------
myMatrix


## ----eval=TRUE,echo=TRUE, fig.height=4,fig.width=6----------------------------
ppm <- myMatrix/colSums(myMatrix)
ppm


## ----eval=TRUE,echo=TRUE, fig.height=4,fig.width=6----------------------------
seqLogo::seqLogo(ppm)


## ----eval=TRUE,echo=TRUE, fig.height=4,fig.width=6----------------------------
Sox9_IC <- toICM(SOX9mat)
TFBSTools::seqLogo(Sox9_IC)


## ----eval=TRUE,echo=TRUE, fig.height=4,fig.width=6----------------------------
library(ggseqlogo)
library(ggplot2)
ggseqlogo(myMatrix)+theme_minimal()


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Motif enrichment analysis

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Motif enrichment analysis

---
"    
  )
  
}


## ----eval=F-------------------------------------------------------------------
# library(Herper)
# install_CondaTools("meme","motif_analysis", channels="bioconda", pathToMiniConda = "miniconda")


## -----------------------------------------------------------------------------

opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["collection"]] <- "CORE"
opts[["all_versions"]] <- FALSE

jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))



## -----------------------------------------------------------------------------

library(universalmotif)
motifs <- getMatrixSet(sq24, opts)
motifs_um <- convert_motifs(motifs)
class(motifs_um)
class(motifs_um[[1]])


## -----------------------------------------------------------------------------

library(rtracklayer)

UpinW6_peaks <- rtracklayer::import("data/UpinW6.bed")



## ----eval=T-------------------------------------------------------------------
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)

peaksSequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, UpinW6_peaks)
names(peaksSequences) <- paste0("peak_", 1:length(UpinW6_peaks))


## ----eval=F, echo=F, warning=F, message=F-------------------------------------
# library(Biostrings)
# mm10 <- readDNAStringSet("~/Downloads/mm10.fa.gz")
# 
# peaksSequences <- getSeq(mm10, UpinW6_peaks)
# names(peaksSequences) <- paste0("peak_", 1:length(UpinW6_peaks))


## ----eval=T-------------------------------------------------------------------
HC_peaks <- rtracklayer::import("data/HC_Peaks.bed")
background_peaksSequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, HC_peaks)
names(background_peaksSequences) <- paste0("peak_", 1:length(HC_peaks))


## ----eval=F, echo=F, warning=F, message=F-------------------------------------
# HC_peaks <- rtracklayer::import("data/HC_Peaks.bed")
# background_peaksSequences <- getSeq(mm10, HC_peaks)
# names(background_peaksSequences) <- paste0("peak_", 1:length(HC_peaks))


## ----eval=F-------------------------------------------------------------------
# library(memes)
# 
# meme_path_id="miniconda/envs/motif_analysis/bin"


## ----eval=F, echo=FALSE-------------------------------------------------------
# 
# meme_path_id ="/Users/mattpaul/Documents/RU/Software/mini/envs/motif_analysis/bin"
# 


## ----eval=F-------------------------------------------------------------------
# ame_results <- runAme(peaksSequences,
#                           control = background_peaksSequences,
#                           outdir = "ame_diff",
#                           meme_path = meme_path_id,
#                           database = motifs_um)


## ----eval=F, echo=F-----------------------------------------------------------
# 
# save(ame_results, file = "data/ame_results.RData")
# 


## -----------------------------------------------------------------------------

load(file = "data/ame_results.RData")

class(ame_results)



## -----------------------------------------------------------------------------

ame_results %>% 
  dplyr::filter(rank %in% 1:10) %>% 
  plot_ame_heatmap()  +
    ggtitle("Top 10 AME Hits")


## -----------------------------------------------------------------------------

top_res_motif <- getMatrixByID(sq24,ame_results$motif_alt_id[1])
top_res_motif_mat <- Matrix(top_res_motif)
ggseqlogo(top_res_motif_mat)+theme_minimal()



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# De novo motif analysis

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# De novo motif analysis

---
"    
  )
  
}


## ----eval=F-------------------------------------------------------------------
# dreme_results <- runDreme(peaksSequences,
#                           control = background_peaksSequences,
#                           outdir = "dreme_diff",
#                           meme_path = meme_path_id)


## ----eval=F, echo=F-----------------------------------------------------------
# 
# save(dreme_results, file = "data/dreme_results.RData")
# 


## -----------------------------------------------------------------------------

load(file = "data/dreme_results.RData")

class(dreme_results)



## -----------------------------------------------------------------------------
dreme_results[1:3,] %>% 
  to_list() %>% 
  view_motifs()


## ----eval=F-------------------------------------------------------------------
# 
# tomtom_results <- runTomTom(dreme_results,
#             meme_path = meme_path_id,
#             database = motifs_um,
#             outdir = "tomtom_from_dreme"
#             )
# 


## ----eval=F, echo=F-----------------------------------------------------------
# 
# save(tomtom_results, file = "data/tomtom_results.RData")
# 


## -----------------------------------------------------------------------------
load(file = "data/tomtom_results.RData")
view_motifs(tomtom_results$best_match_motif[1:3])



## -----------------------------------------------------------------------------
unlist(lapply(tomtom_results$tomtom[[1]]$match_motif, function(x) x@name))



## ----eval=F-------------------------------------------------------------------
# 
# tomtom_results[1,] %>%
#   view_tomtom_hits(top_n = 3)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Finding Motifs

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Finding Motifs

---
"    
  )
  
}


## -----------------------------------------------------------------------------
opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["collection"]] <- "CORE"
opts[["all_versions"]] <- FALSE

jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
motifs <- getMatrixSet(sq24, opts)


## -----------------------------------------------------------------------------
peaksSequences


## -----------------------------------------------------------------------------
library(motifmatchr)
motifs_subset <- motifs[names(motifs) %in% ame_results$motif_alt_id[1:3]] 

motif_positions <- matchMotifs(motifs_subset, peaksSequences[1:100],out="positions")

class(motif_positions)
length(motif_positions)



## -----------------------------------------------------------------------------
motif_positions$MA1152.2


## -----------------------------------------------------------------------------

MA1152.2hits <- motif_positions$MA1152.2
names(MA1152.2hits) <- names(peaksSequences[1:100])
unlist(MA1152.2hits, use.names = TRUE)


## -----------------------------------------------------------------------------
motifHits <- matchMotifs(motifs, peaksSequences, out="matches")
class(motifHits)
motifHits


## -----------------------------------------------------------------------------
mmMatrix <- motifMatches(motifHits)
dim(mmMatrix)
mmMatrix[1:8,1:8]


## ----eval=F, echo=T-----------------------------------------------------------
# totalMotifOccurence <- colSums(mmMatrix)
# totalMotifOccurence[1:4]
# totalMotifOccurence["MA1152.2"]


## ----eval=T, echo=F-----------------------------------------------------------
my_mat<-apply(data.matrix(mmMatrix),2,sum)
names(my_mat)<-colnames(mmMatrix)
my_mat[1:4]

my_mat["MA1152.2"]



## ----eval=T-------------------------------------------------------------------
UpinW6_peaksWithMA1152.2 <- UpinW6_peaks[mmMatrix[,"MA1152.2"] == 1]
UpinW6_peaksWithMA1152.2

