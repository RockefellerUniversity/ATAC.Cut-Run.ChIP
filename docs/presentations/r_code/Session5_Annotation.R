params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T) # delete cache before any merging 



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides != "yes"){
  cat("# CUT&RUN/ATAC (part 5) - Peak annotation and functional enrichment

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Peak Annotation

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Peak Annotation

---
"    
  )
  
}



## ----echo=T, eval=T, echo=T, warning=FALSE,tidy=T,message=FALSE---------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomeInfoDb)
library(ChIPseeker)
library(rtracklayer)

cnrPeaks_GR <- rtracklayer::import("data/peaks/SOX9CNR_W6_rep1_macs_peaks.narrowPeak")


## ----eval=T,echo=T, message=FALSE,messages=FALSE, eval=T, echo=T, warning=FALSE----
peakAnno <- annotatePeak(cnrPeaks_GR, tssRegion=c(-1000, 1000), 
                         TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                         annoDb="org.Mm.eg.db")
class(peakAnno)


## ----eval=T,echo=T, message=F,messages=F, eval=T, echo=T, warning=FALSE,tidy=T----
peakAnno


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------
annotatedPeaksGR <- as.GRanges(peakAnno)
annotatedPeaksDF <- as.data.frame(peakAnno)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------
annotatedPeaksGR[1:2,]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------
annotatedPeaksGR$annotation[1:5]
annotatedPeaksGR$SYMBOL[1:5]


## ----eval=T, echo=T, fig.height=5, fig.width=15, warning=FALSE, tidy=T--------
plotAnnoBar(peakAnno)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,fig.height=5, fig.width=15,tidy=T----
plotDistToTSS(peakAnno)


## ----eval=T, echo=T, fig.height=5, fig.width=15, warning=FALSE, tidy=T--------
library(ggupset)
upsetplot(peakAnno, vennpie=F)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Gene Set Enrichment

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Gene Set Enrichment

---
"    
  )
  
}



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
length(unique(annotatedPeaksGR$geneId))


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
W6MinusD0 <- rio::import("data/W6MinusD0.xlsx")

W6MinusD0[1:5, ]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
W6MinusD0_gr <- GRanges(seqnames = W6MinusD0$seqnames,
                        IRanges(start = W6MinusD0$start, end = W6MinusD0$end),
                        log2FoldChange = W6MinusD0$log2FoldChange,
                        padj = W6MinusD0$padj)

W6MinusD0_gr[1:3, ]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE-----------------------------
W6MinusD0_gr_main <- W6MinusD0_gr[as.vector(seqnames(W6MinusD0_gr)) %in% paste0("chr", c(1:19, "X", "Y", "M"))]
W6MinusD0_gr_up <- W6MinusD0_gr_main[W6MinusD0_gr_main$log2FoldChange > 1 & W6MinusD0_gr_main$padj < 0.05, ]

W6MinusD0_gr_up


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Functional Enrichment with nearest gene

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Functional Enrichment with nearest gene

---
"    
  )
  
}



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T, message = F---------
allGeneGR <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
allGeneGR[1:2,]
allGeneIDs <- allGeneGR$gene_id


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE, message = F, results='hide'----
peakAnno_up <- annotatePeak(W6MinusD0_gr_up, 
                            tssRegion=c(-1000, 1000),
                            TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                            annoDb="org.Mm.eg.db")



## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE, message = F----------------

peakAnno_up


## -----------------------------------------------------------------------------
annotatedPeaksGR_up <- as.GRanges(peakAnno_up)
annotatedPeaksGR_up[1:3]


## -----------------------------------------------------------------------------
annotatedPeaksGR_up_prom <- annotatedPeaksGR_up[grepl("Promoter", annotatedPeaksGR_up$annotation)]
up_promGenes_uniq <- unique(annotatedPeaksGR_up_prom$geneId)

length(up_promGenes_uniq)


## -----------------------------------------------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)

GO_result_prom <- enrichGO(gene = up_promGenes_uniq, 
                      universe = allGeneIDs,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP")
GO_result_prom


## ----eval=T,echo=T, message=F, warning=FALSE,tidy=T---------------------------

up_genes_uniq <- unique(annotatedPeaksGR_up$geneId)

length(up_genes_uniq)


## ----eval=T,echo=T, message=F, warning=FALSE,tidy=T---------------------------

GO_result <- enrichGO(gene = up_genes_uniq, 
                      universe = allGeneIDs,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP")
GO_result


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------
GO_result_df <- data.frame(GO_result, row.names = NULL)
GO_result_df[1:3, ]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T, fig.height=4, fig.width=8----
library(enrichplot)
GO_result_plot <- pairwise_termsim(GO_result)
emapplot(GO_result_plot, showCategory = 20, cex_label_category = 0.7)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------
library(msigdbr)
msigdbr_collections()


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------
library(msigdbr)
msig_t2g_reac <- msigdbr(species = "Mus musculus", 
                    category = "C2", 
                    subcategory = "REACTOME")
msig_t2g_reac <- msig_t2g_reac[ , colnames(msig_t2g_reac) %in% c("gs_name", "entrez_gene")]
msig_t2g_reac[1:3, ]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------
reactome <- enricher(gene = up_genes_uniq, 
                     universe = allGeneIDs,
                     TERM2GENE = msig_t2g_reac)
reactome_df <- data.frame(reactome, row.names = NULL)
# clean up table to print on slide
reactome_df <- reactome_df %>% dplyr::select(-Description, -geneID) %>% dplyr::mutate(ID = substr(ID, 1, 50))
reactome_df[1:5, ]



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Functional Enrichment of distal regions with GREAT

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Functional Enrichment of distal regions with GREAT

---
"    
  )
  
}



## ----eval=T,echo=T, message=F, warning=FALSE,tidy=T---------------------------
library(rGREAT)


## ----eval=T,echo=T, message=F, warning=FALSE,tidy=T---------------------------
great_gobp <- great(W6MinusD0_gr_up, gene_sets = "GO:BP", tss_source = "TxDb.Mmusculus.UCSC.mm10.knownGene")



## ----eval=T,echo=T, message=F, warning=FALSE,tidy=T---------------------------
great_gobp


## ----eval=T,echo=T, message=F, warning=FALSE,tidy=T---------------------------
great_gobp_tab <- getEnrichmentTables(great_gobp)

great_gobp_tab[1:5, ]



## ----eval=T,echo=T, message=F, warning=FALSE,tidy=T---------------------------
# convert to a list of gene sets
reac_gene_sets <- split(msig_t2g_reac$entrez_gene, msig_t2g_reac$gs_name)
reac_gene_sets <- lapply(reac_gene_sets, as.character)  # just to make sure gene IDs are all in character.

reac_gene_sets[1:2]


## ----eval=T,echo=T, message=F, warning=FALSE,tidy=T---------------------------

great_reac <- great(W6MinusD0_gr_up, gene_sets = reac_gene_sets, tss_source = "TxDb.Mmusculus.UCSC.mm10.knownGene", mode = "oneClosest")

great_reac_tab <- getEnrichmentTable(great_reac)
great_reac_tab[1:2, ]


## ----eval=T,echo=T, message=F, warning=FALSE,tidy=T---------------------------
great_genes <- getRegionGeneAssociations(great_reac)

great_genes[1:5]

