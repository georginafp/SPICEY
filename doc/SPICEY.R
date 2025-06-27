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

## ----install, eval=FALSE, echo=TRUE-------------------------------------------
# 
# install.packages("devtools")
# devtools::install_github("georginafp/SPICEY")
# 

## ----show-da-atac, message=FALSE, warning=FALSE, eval=FALSE, echo=TRUE--------
# data("atac")

## ----show-da-rna, message=FALSE, warning=FALSE,eval=FALSE, echo=TRUE----------
# data("rna")

## ----links, message=FALSE, warning=FALSE,eval=FALSE, echo=TRUE----------------
# data("cicero_links")
# head(cicero_links)

## ----retsi, message=FALSE, warning=FALSE, eval=FALSE, echo=TRUE---------------
# retsi <- spicey_retsi(atac)
# head(retsi)

## ----show-retsi, message=FALSE, warning=FALSE, echo = FALSE-------------------
data("retsi")
head(retsi)

## ----getsi, message=FALSE, warning=FALSE, eval=FALSE, echo=TRUE---------------
# getsi <- spicey_getsi(rna)

## ----show-getsi, message=FALSE, warning=FALSE, echo = FALSE-------------------
data("getsi")
head(getsi)

## ----re-gene-nearest, message=FALSE, warning=FALSE, eval=FALSE, echo=TRUE-----
# retsi_gene_nearest <- annotate_with_nearest(retsi)

## ----show-re-gene-nearest, message=FALSE, warning=FALSE, echo=FALSE-----------
data("retsi_gene_nearest")
head(retsi_gene_nearest)

## ----re-gene-coaccessibility, message=FALSE, warning=FALSE, eval=FALSE, echo=TRUE----
# 
# # Filter links for high coaccessibility score
# coacc_links <- cicero_links |>
#   dplyr::filter(coaccess > 0.5)
# 
# # Annotate links with CCANs and filter for Promoter-Distal interactions
# links <- annotate_links_with_ccans(
#   links = coacc_links,
#   coaccess_cutoff_override = 0.25,
#   filter_promoter_distal = TRUE,
#   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
# 
# retsi_gene_coacc <- annotate_with_coaccessibility(
#   re = retsi,
#   links = links,
#   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
#   name_links = "HPAP")
# 

## ----show-re-gene-coacc, message=FALSE, warning=FALSE, echo=FALSE-------------
data("retsi_gene_coacc")
head(retsi_gene_coacc)

## ----spicey-nearest, message=FALSE, warning=FALSE-----------------------------
spicey_nearest <- link_spicey_nearest(retsi_gene_nearest, getsi)
head(spicey_nearest)

## ----spicey-coaccessibility, message=FALSE, warning=FALSE---------------------
spicey_coacc <- link_spicey_coaccessible(retsi_gene_coacc, getsi)
head(spicey_coacc)

## ----spicey-all, message=FALSE, warning=FALSE, eval=FALSE, echo=TRUE----------
# 
# # Compute GETSI
# results <- run_spicey(rna = rna)
# 
# # Compute RETSI
# results <- run_spicey(atac=atac)
# 
# # Compute GETSI + RETSI
# results <- run_spicey(atac=atac, rna=rna)
# 
# # Compute GETSI + RETSI and link RE to target genes through nearest gene method
# results <- run_spicey(rna = rna, atac=atac, link_method = "nearest")
# 
# # Compute GETSI + RETSI and link RE to target genes through coaccessibility method
# results <- run_spicey(
#   atac = atac_data,
#   rna = rna_data,
#   link_method = "coaccessibility",
#   links = coaccessibility_links,
#   coaccess_cutoff_override = 0.25,
#   filter_promoter_distal = TRUE
# )
# 

## ----plot---------------------------------------------------------------------


