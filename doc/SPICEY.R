## ----setup, include = FALSE---------------------------------------------------

library(dplyr)

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

## ----install, eval=FALSE,  echo=TRUE------------------------------------------
# # install.packages("devtools")
# devtools::install_github("georginafp/SPICEY")
# 

## ----scatac, eval=FALSE,  echo=TRUE-------------------------------------------
# 
# merged_df <- readRDS("../data/ATAC_DAR.rds") |> unlist()
# retsi_gr <- regioneR::toGRanges(as.data.frame(merged_df))
# 
# # Annotate with TSS distance
# anno <- ChIPseeker::annotatePeak(retsi_gr,
#   TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
#   verbose = FALSE
# )
# 
# retsi_gr$distanceToTSS <- data.frame(anno)$distanceToTSS
# retsi_gr$annotation <- "Distal"
# retsi_gr$annotation[abs(retsi_gr$distanceToTSS) <= 2000] <- "Promoter"
# 
# # Add nearest gene name
# genes <- biomart_genes()$gr
# retsi_gr$nearestGeneSymbol <- genes$external_gene_name[nearest(retsi_gr, promoters(genes, 1, 0))]
# 
# saveRDS(retsi_gr, "data/FINAL_ATAC.rds")
# 

## ----output-atac, echo=FALSE--------------------------------------------------

head(readRDS("../data/FINAL_ATAC.rds"))


## ----scrna, eval=FALSE, echo=TRUE---------------------------------------------
# 
# gr_list <- readRDS("../data/RNA_DEG.rds")
# 
# gr_df_list <- lapply(gr_list, as.data.frame)
# gr_list <- lapply(names(gr_df_list), function(name) {
#   df <- gr_df_list[[name]]
#   df$symbol <- rownames(df)
#   df$ensembl_id <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
#     keys = df$symbol,
#     column = "ENSEMBL",
#     keytype = "SYMBOL",
#     multiVals = "first"
#   )
# 
#   df_annot <- dplyr::left_join(
#     df,
#     biomart_genes()$df %>% dplyr::select(chromosome_name, start_position, end_position, strand, ensembl_gene_id, gene_biotype),
#     by = c("ensembl_id" = "ensembl_gene_id")
#   )
# 
#   if (nrow(df_annot) == 0) return(NULL)
# 
#   df_annot %>%
#     dplyr::filter(!is.na(chromosome_name),
#                   !is.na(start_position),
#                   !is.na(end_position),
#                   gene_biotype == "protein_coding",
#                   chromosome_name %in% c(as.character(1:22), "X", "Y")) %>%
#     dplyr::mutate(chromosome_name = paste0("chr", chromosome_name)) %>%
#     makeGRangesFromDataFrame(
#       seqnames.field = "chromosome_name",
#       start.field = "start_position",
#       end.field = "end_position",
#       strand.field = "strand",
#       keep.extra.columns = TRUE
#     ) %>%
#     sort()
# })
# 
# gr_list <- Filter(Negate(is.null), gr_list)
# names(gr_list) <- names(gr_df_list)
# 
# # Merge into a single GRanges
# retsi_rna <- purrr::imap(gr_list, ~ {
#   mcols(.x)$cell_type <- .y
#   .x
# }) %>%
#   purrr::reduce(c)
# 
# saveRDS(retsi_rna, "data/FINAL_RNA.rds")
# 

## ----output-rna, echo=FALSE---------------------------------------------------

head(readRDS("../data/FINAL_RNA.rds"))


## ----spicey-nearest, eval=FALSE,  echo=TRUE-----------------------------------
# 
# result_nearest <- run_spicey(
#   atac_path = system.file("extdata", "FINAL_ATAC.rds", package = "SPICEY"),
#   rna_path = system.file("extdata", "FINAL_RNA.rds", package = "SPICEY"),
#   linking_method = "nearest"
# )
# 

## ----spicey-coaccessibility, eval=FALSE,  echo=TRUE---------------------------
# 
# result_coacc <- run_spicey(
#   atac_path = system.file("extdata", "FINAL_ATAC.rds", package = "SPICEY"),
#   rna_path = system.file("extdata", "FINAL_RNA.rds", package = "SPICEY"),
#   links_path = system.file("extdata", "COACC_LINKS.rds", package = "SPICEY"),
#   linking_method = "coaccessibility"
# )
# 

## ----spicey-outputs-----------------------------------------------------------

## ---------------------------------------------------------------------------##
# SPICEY using nearest ---------------------------------------------------------
## ---------------------------------------------------------------------------##
readRDS("../data/SPICEY_nearest.rds") %>% 
  as.data.frame() %>%
  head() %>%
  print()


## ---------------------------------------------------------------------------##
# SPICEY using co-accessibility ------------------------------------------------
## ---------------------------------------------------------------------------##
readRDS("../data/SPICEY_coaccessible.rds") %>% 
  as.data.frame() %>%
  head() %>%
  print()


## ----plot---------------------------------------------------------------------


