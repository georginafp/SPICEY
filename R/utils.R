utils::globalVariables(c(
  "CCAN", "peak", "seqnames", "start", "end", "distanceToTSS", "annotation",
  "annotation1", "annotation2", "avg_FC", "max_FC", "norm_FC", "weight",
  "entropy", "score", "norm_entropy", ".data", "prob", "entropy_component",
  "ccan", "name", "ENTREZID", "GENETYPE", "nearestGeneSymbol", "gene_coacc",
  "gene_id", "GETSI", "cell_type", "GETSI_entropy", "avg_log2FC", "p_val",
  "p_val_adj", "everything"
))


#' Combine Differential Expression Results Across Cell Types
#'
#' This function takes a list of differential expression data frames (from scRNA-seq)
#' for multiple cell types and combines them into a single annotated data frame.
#' Each data frame should contain gene-level statistics (e.g., log2FC, p-value).
#' @param rna_da A named list of data frames. Each element corresponds to a cell type
#' and contains differential expression results with genes as row names.
#' @param gene_id A string indicating the name of the column in each data frame
#' that contains gene identifiers (e.g., gene symbols).
#' @return A single data frame combining all input differential expression results
#' with added columns for gene ids and cell types.
combine_gex_da <- function(rna_da, gene_id) {
  if (is.null(gene_id) || !is.character(gene_id)) {
    stop("You must supply a valid 'gene_id' (character string) to indicate the gene identifier column.")
  }

  gr_list_annot <- lapply(names(rna_da), function(cell_type) {
    df <- rna_da[[cell_type]]
    df$gene_id <- df[[gene_id]]

    if (!"cell_type" %in% colnames(df)) {
      df$cell_type <- cell_type
    }
    df
  })

  names(gr_list_annot) <- names(rna_da)
  gr_list_annot <- Filter(Negate(is.null), gr_list_annot)

  combined_df <- do.call(rbind, gr_list_annot)
  rownames(combined_df) <- NULL
  return(combined_df)
}



#' Get main chromosomes from a GRanges or TxDb, filtering alt/hap/random/etc.
#' @param txdb A TxDb object (e.g., from `GenomicFeatures::makeTxDbFromGFF()`
#' or a prebuilt TxDb package).
#' @param keep_mito Logical, whether to keep mitochondrial/plastid contigs.
#' @param include_only Optional character vector of seqlevels to
#' force-include (manual override).
#' @param verbose Logical, print what is kept/removed.
#' @return Character vector of filtered seqlevels.
get_main_seqlevels <- function(txdb,
                               keep_mito = FALSE,
                               include_only = NULL,
                               verbose = FALSE) {
  seqs <- GenomeInfoDb::seqlevels(txdb)

  if (!is.null(include_only)) {
    return(seqs[seqs %in% include_only])
  }

  exclude_patterns <- c("random", "alt", "fix", "hap", "Un", "decoy", "PATCH", "GL", "KI", "JH", "HSCHR")
  mito_patterns <- c("chrM", "^MT$", "^Mt$", "^Pt$")

  exclude_regex <- paste(exclude_patterns, collapse = "|")
  mito_regex <- paste(mito_patterns, collapse = "|")

  seqs_to_keep <- seqs[!grepl(exclude_regex, seqs, ignore.case = TRUE)]

  if (!keep_mito) {
    seqs_to_keep <- seqs_to_keep[!grepl(mito_regex, seqs_to_keep, ignore.case = TRUE)]
  }

  if (verbose) {
    message("Kept: ", paste(seqs_to_keep, collapse = ", "))
    removed <- setdiff(seqs, seqs_to_keep)
    if (length(removed) > 0) message("Removed: ", paste(removed, collapse = ", "))
  }

  return(seqs_to_keep)
}


#' Get Promoter Regions with Optional Filtering for Protein-Coding Genes
#' Extracts promoter regions (+/- 2kb default) from a TxDb object, optionally restricting to protein-coding genes
#' and retaining only primary chromosomes (using a robust internal filter). Gene symbols and types are added using an AnnotationDbi object.
#' @param txdb A TxDb object (e.g., from `GenomicFeatures::makeTxDbFromGFF()` or a prebuilt TxDb package).
#' @param annot_dbi An AnnotationDbi object (e.g., `org.Hs.eg.db`, `org.Mm.eg.db`) for mapping gene IDs to symbols and types.
#' @param keep_mito Logical; whether to keep mitochondrial chromosome (default = FALSE).
#' @param protein_coding_only Logical; whether to restrict to protein-coding genes (default = TRUE).
#' @param verbose Logical; print informative messages (default = TRUE).
#' @return A `GRanges` object of promoter regions, annotated with gene symbols, optionally filtered to protein-coding genes.
#' @importFrom GenomicFeatures transcripts promoters genes
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom dplyr filter pull
get_promoters <- function(txdb,
                          annot_dbi,
                          keep_mito = FALSE,
                          protein_coding_only = TRUE,
                          verbose = TRUE) {
  # Get main chromosomes using helper (no unused args here)
  main_chrs <- get_main_seqlevels(txdb,
                                  keep_mito = keep_mito,
                                  verbose = verbose)

  # Fetch promoters
  if (length(main_chrs) < 3) {
    if (verbose) message("Few chromosomes after filtering (", length(main_chrs), ") - skipping keepSeqlevels()")
    proms <- GenomicFeatures::promoters(GenomicFeatures::transcripts(txdb, columns = "GENEID"),
                                        upstream = 2000, downstream = 2000)
  } else {
    proms <- GenomicFeatures::transcripts(txdb, columns = "GENEID") |>
      GenomicRanges::promoters(upstream = 2000, downstream = 2000) |>
      GenomeInfoDb::keepSeqlevels(main_chrs, pruning.mode = "coarse")
  }

  # Replace empty GENEID with NA
  proms$GENEID <- vapply(proms$GENEID, function(x) if (length(x) == 0) NA_character_ else x, FUN.VALUE = character(1))


  # Add gene symbols
  proms$gene_id <- AnnotationDbi::select(annot_dbi,
                                         keys = proms$GENEID,
                                         columns = "SYMBOL",
                                         keytype = "ENTREZID",
                                         multiVals = "first")$SYMBOL

  # Optionally filter to protein-coding genes
  if (protein_coding_only) {
    pc <- AnnotationDbi::select(annot_dbi,
                                keys = proms$GENEID,
                                columns = "GENETYPE",
                                keytype = "ENTREZID",
                                multiVals = "first") |>
      dplyr::filter(!is.na(ENTREZID), !is.na(GENETYPE), GENETYPE == "protein-coding") |>
      dplyr::pull(ENTREZID)

    proms <- proms[which(proms$GENEID %in% pc)]
  }

  return(proms)
}
