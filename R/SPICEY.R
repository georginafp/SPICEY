#' SPICEY for tissue specificity analysis
#'
#' Computes tissue-specificity scores from differential accessibility (RETSI)
#' and/or gene expression (GETSI) data obtained from single cell experiments. Supports:
#' \itemize{
#'   \item RETSI calculation from differential accessibility data in different cell types/clusters (scATAC-seq).
#'   \item GETSI calculation from differential expression data in different cell types/clusters (scRNA-seq).
#'   \item Optional integration of RETSI and GETSI scores by linking gene associations
#'   (see \code{\link{annotate_with_nearest}} or \code{\link{annotate_with_coaccessibility}})
#'   }
#' @param rna Either a single \code{data.frame} or a named list of
#'   \code{data.frame}s or \code{GRanges} where each element corresponds
#'   to a cell type. It should contain differential expression results,
#'   with required columns:
#'   \describe{
#'     \item{gene_id}{Identifier of the gene. Must be an official gene symbol (e.g., \code{GAPDH}).
#'     The name of this column should match the \code{gene_id} argument.}
#'     \item{avg_log2FC}{Average log2 fold-change for the gene in that cell type.}
#'     \item{p_val_adj}{Adjusted p-value (e.g., FDR-corrected).}
#'     \item{cell_type}{Cell type or cluster label. Only necessary when input is
#'           a single \code{data.frame}. If input is a list, it will be generated from list names
#'   Note that the same gene may appear multiple times across cell types.}
#'   }
#' @param atac Either a single \code{data.frame} or a named list of
#'   \code{data.frame}s or \code{GRanges} where each element corresponds
#'   to a cell type.  It should contain differential chromatin accessibility
#'   results with required columns:
#'   \describe{
#'     \item{region_id}{Unique identifier of the region (e.g., \code{chr1-5000-5800}).}
#'     \item{avg_log2FC}{Average log2 fold-change for accessibility in that cell type.}
#'     \item{p_val_adj}{Adjusted p-value (e.g., FDR-corrected).}
#'     \item{cell_type}{Cell type or cluster label. Only necessary when input is
#'           a single \code{data.frame}. If input is a list, it will be generated
#'           from list names.
#'   Note that the same region may appear multiple times across cell types.}
#'   }
#' @param annotation (Optional). A data.frame linking \code{gene_id} to \code{region_id}.
#' They should have the same names provided in the respective parameters.
#' This can be provided by the user or generated using the functions:
#' \code{\link{annotate_with_nearest}} or \code{\link{annotate_with_coaccessibility}}.
#' It should contain at least the following columns:
#' \describe{
#'   \item{region_id}{Unique identifier of the region (e.g., \code{chr1-5000-5800}).
#'     The name of this column should match the \code{region_id} argument.}
#'   \item{cell_type}{Cell type or cluster label. (e.g., \code{Acinar})}
#'   \item{gene_id}{Identifier of the gene. Must be an official gene symbol (e.g., \code{GAPDH}).
#'     The name of this column should match the \code{gene_id} argument.}}
#' @param verbose Logical; print messages (default TRUE).
#' @return Depending on inputs, returns RETSI and/or GETSI data frames, optionally linked and annotated.
#' @examples
#' data(rna)
#' data(atac)
#'
#' # Calculate RETSI only
#' retsi <- SPICEY(atac = atac)
#'
#' # Calculate GETSI only
#' getsi <- SPICEY(rna = rna)
#'
#' # Calculate both
#' both <- SPICEY(
#'   rna = rna,
#'   atac = atac
#' )
#'
#' # Integrate RETSI and GETSI with nearest gene
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#' peaks <- unique(unlist(atac)[, c("region_id")])
#'
#' annotation_near <- annotate_with_nearest(
#'   peaks = peaks,
#'   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   annot_dbi = org.Hs.eg.db,
#'   protein_coding_only = TRUE,
#'   verbose = TRUE,
#'   add_tss_annotation = FALSE,
#'   upstream = 2000,
#'   downstream = 2000
#' )
#' spicey_near <- SPICEY(
#'   rna = rna,
#'   atac = atac,
#'   annotation = annotation_near
#' )
#'
#' # Integrate RETSI and GETSI with coaccessibility
#' data(cicero_links)
#'
#' annotation_coacc <- annotate_with_coaccessibility(
#'   peaks = peaks,
#'   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   links_df = cicero_links,
#'   annot_dbi = org.Hs.eg.db,
#'   protein_coding_only = TRUE,
#'   verbose = TRUE,
#'   add_tss_annotation = FALSE,
#'   upstream = 2000,
#'   downstream = 2000
#' )
#' spicey_coacc <- SPICEY(
#'   rna = rna,
#'   atac = atac,
#'   annotation = annotation_coacc
#' )
#' @export
SPICEY <- function(atac = NULL,
                   rna = NULL,
                   annotation = NULL,
                   verbose = TRUE) {
  if (is.null(atac) && is.null(rna)) {
    stop("Provide at least one of 'atac' or 'rna'.")
  }

  if (!is.null(rna)) {
    if (verbose) message("Computing GETSI & entropy...")
    rna <- .parse_input_diff(rna)

    if (!("gene_id" %in% colnames(rna))) stop("'gene_id' column should be present in rna.")

    getsi <- compute_spicey_index(diff = rna, id = "gene_id") |>
      dplyr::rename(GETSI = score, GETSI_entropy = norm_entropy)
  } else {
    getsi <- NULL
  }
  if (!is.null(atac)) {
    if (verbose) message("Computing RETSI & entropy...")
    atac <- .parse_input_diff(atac)

    if (!("region_id" %in% colnames(atac))) stop("'region_id' column should be present in atac.")


    retsi <- compute_spicey_index(atac, id = "region_id") |>
      dplyr::rename(RETSI = score, RETSI_entropy = norm_entropy)
  } else {
    retsi <- NULL
  }
  results <- list(RETSI = retsi, GETSI = getsi)
  results <- results[!vapply(results, is.null, logical(1))]

  if (is.null(annotation) & length(results) > 1) {
    if (verbose) message("SPICEY pipeline successfully completed")
    return(results)
  } else if (is.null(annotation) & length(results) == 1) {
    if (verbose) message("SPICEY pipeline successfully completed")
    return(results[[1]])
  } else {
    if (verbose) message("Linking RETSI and GETSI using provided annotation...")
    combined <- link_spicey(
      retsi = retsi, getsi = getsi, annotation = annotation
    )
    results$linked <- combined
    if (verbose) message("SPICEY pipeline successfully completed")
    return(results)
  }
}
