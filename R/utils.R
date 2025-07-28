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
#' Merges differential expression results from multiple cell types into one data frame.
#' Designed for scRNA-seq outputs where each cell type is analyzed separately.
#' @param rna List of data frames (or GRanges-like objects), each with differential expression results per cell type.
#'   Each must include:
#'   \describe{
#'     \item{\code{[gene_id]}}{Gene identifier column matching \code{gene_id} argument.}
#'     \item{\code{avg_log2FC}}{Average log2 fold-change.}
#'     \item{\code{p_val}}{Raw p-value.}
#'     \item{\code{p_val_adj}}{Adjusted p-value.}
#'   }
#'   Adds \code{cell_type} column if missing, using list names.
#' @param gene_id Column name for gene identifiers.
#' @return A single \code{data.frame} combining all results with columns:
#' \itemize{
#'   \item \code{gene_id}: Renamed from the specified gene identifier column.
#'   \item \code{cell_type}: Derived from the list element names (or original column if present).
#'   \item All original differential expression metrics (\code{avg_log2FC}, \code{p_val},
#'   \code{p_val_adj}).
#' }
combine_gex_da <- function(rna, gene_id) {
  if (is.null(gene_id) || !is.character(gene_id)) {
    stop("You must supply a valid 'gene_id' (character string) to indicate the gene identifier column.")
  }

  gr_list_annot <- lapply(names(rna), function(cell_type) {
    df <- rna[[cell_type]]
    df$gene_id <- df[[gene_id]]

    if (!"cell_type" %in% colnames(df)) {
      df$cell_type <- cell_type
    }
    df
  })

  names(gr_list_annot) <- names(rna)
  gr_list_annot <- Filter(Negate(is.null), gr_list_annot)
  combined_df <- do.call(rbind, gr_list_annot)
  rownames(combined_df) <- NULL
  return(combined_df)
}



#' Filter Main Chromosomes from Genome Annotation
#'
#' Retrieves primary chromosomes from a genome annotation, excluding alternative contigs,
#' haplotypes, decoys, and patches. Optionally keeps mitochondrial chromosomes or includes
#' specific seqlevels.
#' @param txdb A \code{TxDb}, \code{GRanges}, or any object supporting \code{seqlevels()}.
#' @param keep_mito Logical (default = \code{FALSE}). Keep mitochondrial chromosomes.
#' @param include_only Character vector of seqlevels to keep (ignores other filters if set).
#' @param verbose Logical (default = \code{FALSE}). Print retained/removed seqlevels.
#' @return Character vector of filtered seqlevels.
#' @details Excludes seqlevels matching: \code{"random"}, \code{"alt"}, \code{"fix"}, \code{"hap"},
#' \code{"Un"}, \code{"decoy"}, \code{"PATCH"}, \code{"GL"}, \code{"KI"}, \code{"JH"}, \code{"HSCHR"}.
#' Mitochondrial chromosomes match: \code{"chrM"}, \code{"^MT$"}, \code{"^Mt$"}, \code{"^Pt$"}.
get_main_seqlevels <- function(txdb,
                               keep_mito = FALSE,
                               include_only = NULL,
                               verbose = FALSE) {
  seqs <- GenomeInfoDb::seqlevels(txdb)

  if (!is.null(include_only)) {
    return(seqs[seqs %in% include_only])
  }

  exclude_patterns <- c("random", "alt", "fix", "hap", "Un", "decoy", "PATCH",
                        "GL", "KI", "JH", "HSCHR")
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




#' Extract Promoter Regions from TxDb
#'
#' Extracts promoter regions (default: +/-2 kb from TSS) from a \code{TxDb} object,
#' optionally filtering for protein-coding genes and main chromosomes. Gene symbols
#' and biotypes are added using an \code{AnnotationDbi} object.
#' @param txdb \code{TxDb} object from \pkg{GenomicFeatures}.
#' @param annot_dbi \code{AnnotationDbi} object (e.g., \pkg{org.Hs.eg.db}) for gene symbols/biotypes.
#' @param keep_mito Logical (default = \code{FALSE}). Keep mitochondrial chromosomes.
#' @param protein_coding_only Logical (default = \code{TRUE}). Restrict to protein-coding genes.
#' @param verbose Logical (default = \code{TRUE}). Print filtering details.
#' @return A \code{GRanges} of promoters with \code{gene_id}, \code{symbol}, and \code{gene_type}.
#' @importFrom GenomicFeatures transcripts promoters genes
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom dplyr filter pull
get_promoters <- function(txdb,
                          annot_dbi,
                          keep_mito = FALSE,
                          protein_coding_only = TRUE,
                          verbose = TRUE) {

  main_chrs <- get_main_seqlevels(txdb,
                                  keep_mito = keep_mito,
                                  verbose = verbose)

  if (length(main_chrs) < 3) {
    if (verbose) message("Few chromosomes after filtering (", length(main_chrs), ") - skipping keepSeqlevels()")
    proms <- GenomicFeatures::promoters(GenomicFeatures::transcripts(txdb, columns = "GENEID"),
                                        upstream = 2000, downstream = 2000)
  } else {
    proms <- GenomicFeatures::transcripts(txdb, columns = "GENEID") |>
      GenomicRanges::promoters(upstream = 2000, downstream = 2000) |>
      GenomeInfoDb::keepSeqlevels(main_chrs, pruning.mode = "coarse")
  }

  proms$GENEID <- vapply(proms$GENEID, function(x) if (length(x) == 0) NA_character_ else x, FUN.VALUE = character(1))
  proms$gene_id <- AnnotationDbi::select(annot_dbi,
                                         keys = proms$GENEID,
                                         columns = "SYMBOL",
                                         keytype = "ENTREZID",
                                         multiVals = "first")$SYMBOL

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
