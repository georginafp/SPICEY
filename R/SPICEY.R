#' Runs the full SPICEY pipeline
#'
#' This function runs the full SPICEY pipeline to compute tissue specificity
#' scores from single cell ATAC-seq (RETSI) and/or RNA-seq (GETSI) data.
#' It annotates regulatory elements to target genes via nearest-gene or
#' co-accessibility annotation, and optionally integrates RETSI and GETSI measures
#' by linking them using the specified annotation method.
#'
#' @param atac Single-cell ATAC-seq data (Seurat or SummarizedExperiment).
#' @param rna Single-cell RNA-seq data (Seurat or SummarizedExperiment).
#' @param gene_id Column name for gene IDs; required if RNA data is provided.
#' @param annot_method Character; Annotation method to link regulatory elements to genes.
#'   Must be either `"nearest"` or `"coaccessibility"`. This specifies which
#'   gene annotation strategy is used.
#' @param links Co-accessibility links data; required if `annot_method` is `"coaccessibility"`.
#' @param link_spicey_measures Logical; whether to link RETSI and GETSI by gene using `link_spicey`.
#' @param coaccess_cutoff_override Numeric cutoff for co-accessibility (default 0.25).
#' @param filter_promoter_distal Logical; filter promoter-distal links (default TRUE).
#' @param filter_protein_coding Logical; filter to protein-coding genes (default TRUE).
#' @param keep_mito Logical; whether to keep mitochondrial regions (default FALSE).
#' @param txdb TxDb object for annotation; required if annotation is used.
#' @param annot_dbi AnnotationDbi object for gene info; required if annotation is used.
#' @param add_tss_annotation Logical; include TSS annotation columns (default FALSE).
#' @param verbose Logical; print messages during execution (default TRUE).
#' @return Either a tibble/data.frame or a list:
#'   If \code{link_spicey_measures} is \code{TRUE}, returns a combined tibble/data.frame
#'   with specificity and annotation values. If \code{FALSE}, returns a list containing
#'   separate elements for RETSI and GETSI, each with their respective data.
#' @export
#' @examples
#' data(atac)
#' data(rna)
#' data(cicero_links)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' spicey_nearest <- SPICEY(
#'   atac = atac,
#'   rna = rna,
#'   gene_id = "gene_id",
#'   annot_method = "nearest",
#'   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
#'   link_spicey_measures = TRUE
#' )
#'
#' spicey_coacc <- SPICEY(
#'   atac = atac,
#'   rna = rna,
#'   gene_id = "gene_id",
#'   annot_method = "coaccessibility",
#'   links = cicero_links,
#'   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
#'   link_spicey_measures = TRUE,
#'   coaccess_cutoff_override = 0.25
#' )
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

    # Use the generalized link_spicey function with explicit method
    combined <- link_spicey(
      retsi_annotated = retsi_annotated,
      getsi = getsi,
      method = annot_method
    )

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
