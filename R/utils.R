#' Get main chromosomes from a GRanges or TxDb, filtering alt/hap/random/etc.
#'
#' @param x A TxDb, GRanges, or similar object.
#' @param keep_mito Logical, whether to keep mitochondrial/plastid contigs.
#' @param include_only Optional character vector of seqlevels to force-include (manual override).
#' @param verbose Logical, print what is kept/removed.
#' @return Character vector of filtered seqlevels.
#' @export
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
    message("→ Kept: ", paste(seqs_to_keep, collapse = ", "))
    removed <- setdiff(seqs, seqs_to_keep)
    if (length(removed) > 0) message("→ Removed: ", paste(removed, collapse = ", "))
  }

  return(seqs_to_keep)
}


#' Get Promoter Regions with Optional Filtering for Protein-Coding Genes
#' Extracts promoter regions (±2kb default) from a TxDb object, optionally restricting to protein-coding genes
#' and retaining only primary chromosomes (using a robust internal filter). Gene symbols and types are added using an AnnotationDbi object.
#' @param txdb A TxDb object (e.g., from `GenomicFeatures::makeTxDbFromGFF()` or a prebuilt TxDb package).
#' @param annot_dbi An AnnotationDbi object (e.g., `org.Hs.eg.db`, `org.Mm.eg.db`) for mapping gene IDs to symbols and types.
#' @param keep_mito Logical; whether to keep mitochondrial chromosome (default = FALSE).
#' @param protein_coding_only Logical; whether to restrict to protein-coding genes (default = TRUE).
#' @param verbose Logical; print informative messages (default = TRUE).
#' @return A `GRanges` object of promoter regions, annotated with gene symbols, optionally filtered to protein-coding genes.
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import GenomeInfoDb
#' @importFrom AnnotationDbi select
#' @importFrom dplyr filter pull
#' @export
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
    if (verbose) message("Few chromosomes after filtering (", length(main_chrs), ") — skipping keepSeqlevels()")
    proms <- GenomicFeatures::promoters(GenomicFeatures::transcripts(txdb, columns = "GENEID"),
                                        upstream = 2000, downstream = 2000)
  } else {
    proms <- GenomicFeatures::transcripts(txdb, columns = "GENEID") |>
      GenomicRanges::promoters(upstream = 2000, downstream = 2000) |>
      GenomeInfoDb::keepSeqlevels(main_chrs, pruning.mode = "coarse")
  }

  # Replace empty GENEID with NA
  proms$GENEID <- sapply(proms$GENEID, function(x) ifelse(length(x) == 0, NA, x))

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
