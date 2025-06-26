## ----setup, include = FALSE---------------------------------------------------

library(dplyr)
library(GenomicRanges)
library(cicero)

knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    eval = TRUE,
    warning = FALSE,
    message = FALSE,
    fig.align = "center",
    out.width = "80%"
)

library(SPICEY)

## ----logo, echo=FALSE, eval=TRUE, out.width='10%'-----------------------------
knitr::include_graphics("../man/figures/logo_spicey.png", dpi = 800)

## ----install, eval=FALSE,  echo=TRUE------------------------------------------
# 
# install.packages("devtools")
# devtools::install_github("georginafp/SPICEY")
# 

## ----show-da-atac, message=FALSE, warning=FALSE-------------------------------
data("atac")

## ----show-da-rna, message=FALSE, warning=FALSE--------------------------------
data("rna")

## ----show-links, message=FALSE, warning=FALSE---------------------------------
data("cicero_links")

## ----retsi, message=FALSE, warning=FALSE--------------------------------------
retsi <- spicey_retsi(atac)
head(retsi)

## ----getsi, message=FALSE, warning=FALSE--------------------------------------
getsi <- spicey_getsi(rna)
head(getsi)

## ----re-gene-nearest, message=FALSE, warning=FALSE----------------------------
retsi_gene_nearest <- annotate_with_nearest(retsi)
head(retsi_gene_nearest)

## ----re-gene-coaccessibility, message=FALSE, warning=FALSE--------------------

coacc_links <- cicero_links |> 
  dplyr::filter(coaccess > 0.5)

links <- annotate_links_with_ccans(
  links = coacc_links,
  coaccess_cutoff_override = 0.25,
  filter_promoter_distal = TRUE,
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)

retsi_gene_coacc <- annotate_with_coaccessibility(
  re = retsi,
  links = links,
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
  name_links = "HPAP")

head(retsi_gene_coacc)

## ----spicey-nearest, message=FALSE, warning=FALSE-----------------------------
spicey_nearest <- link_spicey_nearest(retsi_gene_nearest, getsi)
head(spicey_nearest)

## ----spicey-coaccessibility, message=FALSE, warning=FALSE---------------------
spicey_coacc <- link_spicey_coaccessible(retsi_gene_coacc, getsi)
head(spicey_coacc)

## ----plot---------------------------------------------------------------------


