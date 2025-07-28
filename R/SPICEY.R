#' Run the full SPICEY pipeline for tissue specificity analysis
#'
#' Computes tissue-specificity scores from single-cell chromatin accessibility (RETSI)
#' and/or gene expression (GETSI) data. Supports:
#' \itemize{
#'   \item RETSI calculation from scATAC-seq data.
#'   \item GETSI calculation from scRNA-seq data.
#'   \item Optional annotation of regulatory elements to genes via \code{"nearest"} or \code{"coaccessibility"} methods.
#'   \item Optional integration of RETSI and GETSI scores by linking gene associations.
#' }
#' @param rna A list of data frames (or \code{GRanges}-like objects).
#'   Each element corresponds to a cell type containing differential expression
#'   results, with required columns:
#'   \describe{
#'     \item{gene_id}{Identifier of the gene (e.g., gene symbol, Ensembl ID).}
#'     \item{avg_log2FC}{Average log2 fold-change for the gene in that cell type.}
#'     \item{p_val}{Raw p-value for the differential test.}
#'     \item{p_val_adj}{Adjusted p-value (e.g., FDR-corrected).}
#'     \item{cell_type}{Cell type or cluster label.}
#'   }
#'   Note that the same peak may appear multiple times across cell types.
#' @param atac A list of data frames (or \code{GRanges}-like objects).
#'   Each element corresponds to a cell type containing differential chromatin
#'   accessibility results with required columns:
#'   \describe{
#'     \item{seqnames}{Chromosome name of the regulatory region (e.g., "chr1").}
#'     \item{start}{Start coordinate of the peak.}
#'     \item{end}{End coordinate of the peak.}
#'     \item{avg_log2FC}{Average log2 fold-change for accessibility in that cell type.}
#'     \item{p_val}{Raw p-value for the differential test.}
#'     \item{p_val_adj}{Adjusted p-value (e.g., FDR).}
#'     \item{cell_type}{Cell type or cluster name.}
#'   }
#' @param gene_id A character string specifying the column name in each list element
#'   that contains the gene identifiers.
#' @param annot_method Optional annotation method: \code{"nearest"} or \code{"coaccessibility"}.
#' @param links \code{GInteractions} object with co-accessibility links (required if \code{annot_method = "coaccessibility"}).
#' @param link_spicey_measures Logical; link RETSI and GETSI scores (default FALSE).
#' @param coaccess_cutoff_override Numeric; cutoff for co-accessibility clustering (default 0.25).
#' @param filter_promoter_distal Logical; filter to promoter-distal links (default TRUE).
#' @param filter_protein_coding Logical; restrict to protein-coding genes (default TRUE).
#' @param keep_mito Logical; keep mitochondrial chromosomes (default FALSE).
#' @param txdb \code{TxDb} object for genome annotation (required if annotation requested).
#' @param annot_dbi \code{AnnotationDbi} object for gene ID mapping (required if annotation requested).
#' @param add_tss_annotation Logical; annotate regulatory elements overlapping TSS (default FALSE).
#' @param verbose Logical; print messages (default TRUE).
#' @return Depending on inputs, returns RETSI and/or GETSI data frames, optionally linked and annotated.
#' @examples
#' data(atac); data(rna); data(cicero_links)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'
#' SPICEY(atac, rna, gene_id = "gene_id", annot_method = "nearest",
#'        txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'        annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
#'        link_spicey_measures = TRUE)
#'
#' SPICEY(atac, rna, gene_id = "gene_id", annot_method = "coaccessibility",
#'        links = cicero_links,
#'        txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'        annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
#'        link_spicey_measures = TRUE,
#'        coaccess_cutoff_override = 0.25)
#' @importFrom GenomicRanges promoters distanceToNearest makeGRangesFromDataFrame mcols
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels
#' @importFrom tidyr unnest
#' @importFrom dplyr select
#' @importFrom S4Vectors queryHits subjectHits
#' @export
SPICEY <- function(atac = NULL,
                   rna = NULL,
                   gene_id = NULL,
                   annot_method = NULL,
                   links = NULL,
                   link_spicey_measures = FALSE,
                   coaccess_cutoff_override = 0.25,
                   filter_promoter_distal = TRUE,
                   filter_protein_coding = TRUE,
                   txdb = NULL,
                   keep_mito = FALSE,
                   annot_dbi = NULL,
                   add_tss_annotation = FALSE,
                   verbose = TRUE) {

  if (is.null(atac) && is.null(rna)) {
    stop("Provide at least one of 'atac' or 'rna'.")
  }

  if (!is.null(rna) && is.null(gene_id)) {
    stop("'gene_id' is required when RNA data is supplied.")
  }

  results <- compute_spicey_index(atac = atac, rna = rna, gene_id = gene_id)

  if (is.list(results) && all(c("retsi", "getsi") %in% names(results))) {
    retsi <- results$retsi
    getsi <- results$getsi
  } else {
    if (!is.null(atac) && is.null(rna)) {
      retsi <- results
      getsi <- NULL
    } else if (is.null(atac) && !is.null(rna)) {
      retsi <- NULL
      getsi <- results
    } else {
      retsi <- NULL
      getsi <- NULL
    }
  }

  if (is.null(annot_method)) {
    if (!is.null(retsi) && is.null(getsi)) return(retsi)
    if (is.null(retsi) && !is.null(getsi)) return(getsi)
    return(results)
  }

  if (is.null(txdb)) stop("Annotation requires a 'txdb' object.")
  if (is.null(annot_dbi)) stop("Annotation requires an 'annot_dbi' object.")
  if (is.null(retsi)) stop("ATAC data (RETSI) is required for region-to-gene annotation when 'annot_method' is specified.")

  retsi_annotated <- switch(
    annot_method,
    nearest = annotate_with_nearest(
      retsi = retsi,
      txdb = txdb,
      keep_mito = keep_mito,
      annot_dbi = annot_dbi,
      protein_coding_only = filter_protein_coding,
      verbose = verbose,
      add_tss_annotation = add_tss_annotation
    ),
    coaccessibility = {
      if (is.null(links)) stop("'links' must be provided for coaccessibility annotation.")
      annotate_with_coaccessibility(
        links = links,
        retsi = retsi,
        txdb = txdb,
        keep_mito = keep_mito,
        annot_dbi = annot_dbi,
        protein_coding_only = filter_protein_coding,
        verbose = verbose,
        coaccess_cutoff_override = coaccess_cutoff_override,
        filter_promoter_distal = filter_promoter_distal,
        add_tss_annotation = add_tss_annotation
      )
    },
    stop("Invalid 'annot_method'. Choose 'nearest' or 'coaccessibility'.")
  )

  if (link_spicey_measures) {
    if (is.null(getsi)) stop("RNA data must be provided to link RETSI and GETSI.")
    combined <- link_spicey(
      retsi_annotated = retsi_annotated,
      getsi = getsi,
      method = annot_method)
    message("SPICEY pipeline successfully completed")
    return(combined)
  }

  out <- list(retsi_annotated = retsi_annotated)
  if (!is.null(getsi)) out$getsi <- getsi
  message("SPICEY pipeline successfully completed")

  if (length(out) == 1) {
    return(out[[1]])
  } else {
    return(out)
  }
}
