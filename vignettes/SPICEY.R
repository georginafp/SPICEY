## ----setup, include = FALSE---------------------------------------------------

library(dplyr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(SPICEY)

knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    eval = TRUE,
    warning = FALSE,
    message = FALSE,
    fig.align = "center",
    out.width = "80%"
)


## ----logo, echo=FALSE, eval=TRUE, out.width='10%'-----------------------------
knitr::include_graphics("../man/figures/logo_spicey.png", dpi = 800)

## ----install, eval=FALSE, echo=TRUE-------------------------------------------
# install.packages("devtools")
# devtools::install_github("georginafp/SPICEY")

## ----da-atac, message=FALSE, warning=FALSE------------------------------------
data("atac")

## ----da-rna, message=FALSE, warning=FALSE-------------------------------------
data("rna")

## ----links, message=FALSE, warning=FALSE--------------------------------------
data("cicero_links")
head(cicero_links)

## ----retsi, message=FALSE, warning=FALSE--------------------------------------
retsi <- SPICEY(atac=atac, region_id = "region_id")
head(retsi)

## ----getsi, message=FALSE, warning=FALSE--------------------------------------
getsi <- SPICEY(rna=rna, gene_id = "gene_id")
head(getsi)

## ----spicey, message=FALSE, warning=FALSE-------------------------------------
spicey <- SPICEY(atac=atac, 
                region_id = "region_id", 
                rna=rna,
                gene_id = "gene_id")
lapply(spicey, head)

## ----re-gene-nearest, message=FALSE, warning=FALSE----------------------------

peaks <- SPICEY:::.parse_input_diff(atac)
peaks <- peaks %>% tidyr::separate(region_id,
                into = c("chr", "start", "end"), sep = "-",
                convert = TRUE,remove = FALSE) %>%
                GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

annotation_near <- annotate_with_nearest(peaks = peaks,
                                         txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                         annot_dbi = org.Hs.eg.db,
                                         protein_coding_only = TRUE,
                                         verbose = TRUE,
                                         add_tss_annotation = FALSE,
                                         upstream = 2000,
                                         downstream = 2000)
head(annotation_near)

## ----re-gene-coaccessibility, message=FALSE, warning=FALSE--------------------

peaks <- SPICEY:::.parse_input_diff(atac)
peaks <- peaks %>% tidyr::separate(region_id,
                into = c("chr", "start", "end"), sep = "-",
                convert = TRUE,remove = FALSE) %>%
                GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

annotation_coacc <- annotate_with_coaccessibility(peaks = peaks,
                                                  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                                  links_df=cicero_links,
                                                  annot_dbi = org.Hs.eg.db,
                                                  protein_coding_only = TRUE,
                                                  verbose = TRUE,
                                                  add_tss_annotation = FALSE,
                                                  upstream = 2000,
                                                  downstream = 2000)
head(annotation_coacc)

