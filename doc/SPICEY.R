## ----setup, include = FALSE---------------------------------------------------

library(dplyr)
library(GenomicRanges)
library(cicero)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


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
retsi <- spicey_retsi(atac)
head(retsi)

## ----getsi, message=FALSE, warning=FALSE--------------------------------------
getsi <- spicey_getsi(rna)
head(getsi)

## ----re-gene-nearest, message=FALSE, warning=FALSE----------------------------
retsi_gene_nearest <- annotate_with_nearest(
  retsi = retsi,
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
  annot_dbi = org.Hs.eg.db::org.Hs.eg.db
)

head(retsi_gene_nearest)

## ----re-gene-coaccessibility, message=FALSE, warning=FALSE--------------------

# Filter links for high coaccessibility score
coacc_links <- cicero_links |> 
  dplyr::filter(coaccess > 0.5)

retsi_gene_coacc <- annotate_with_coaccessibility(
  links = coacc_links,
  retsi = retsi,
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
  annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
  coaccess_cutoff_override = 0.25
  )

head(retsi_gene_coacc)

## ----spicey-nearest, message=FALSE, warning=FALSE-----------------------------
spicey_nearest <- link_spicey_nearest(retsi_gene_nearest, getsi)
head(spicey_nearest)

## ----spicey-coaccessibility, message=FALSE, warning=FALSE---------------------
spicey_coacc <- link_spicey_coaccessible(retsi_gene_coacc, getsi)
head(spicey_coacc)

## ----spicey-all, message=FALSE, warning=FALSE---------------------------------

spicey_nearest <- run_spicey(
  atac = atac, 
  rna = rna, 
  annot_method = "nearest", 
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
  annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
  link_spicey_measures = TRUE
)
head(spicey_nearest)

spicey_coacc <- run_spicey(
  atac = atac,
  rna = rna,
  annot_method = "coaccessibility",
  links = coacc_links,
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
  annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
  coaccess_cutoff_override = 0.25,
  link_spicey_measures = TRUE
)
head(spicey_coacc)

## ----plot---------------------------------------------------------------------


