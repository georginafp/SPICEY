#' Run the full SPICEY pipeline
#'
#' This function runs the SPICEY pipeline to compute tissue specificity scores
#' from single-cell ATAC-seq (RETSI) and/or RNA-seq (GETSI) data. It optionally
#' links regulatory elements to genes using either nearest-gene or co-accessibility
#' based annotation, and can combine RETSI and GETSI values for downstream integration.
#'
#' @param atac A Seurat or SummarizedExperiment object containing scATAC-seq data.
#' @param rna A Seurat or SummarizedExperiment object containing scRNA-seq data.
#' @param gene_id A string indicating the name of the column in each data frame that contains gene identifiers
#'   (e.g., gene symbols). Required if \code{rna} is provided.
#' @param annot_method Method for region-to-gene annotation: "nearest" or "coaccessibility".
#' @param links Co-accessibility links (e.g., from Cicero or ArchR) for co-accessibility annotation.
#' @param link_spicey_measures Logical; whether to link RETSI and GETSI via annotated genes.
#' @param coaccess_cutoff_override Minimum co-accessibility score to retain links (default: 0.25).
#' @param filter_promoter_distal Logical; whether to exclude promoter-proximal links (default: TRUE).
#' @param filter_protein_coding Logical; whether to filter to protein-coding genes only (default: TRUE).
#' @param txdb A TxDb object for genomic annotations.
#' @param keep_mito Logical; whether to keep mitochondrial regions (default: FALSE).
#' @param annot_dbi Optional AnnotationDbi object for additional gene info (e.g., gene symbols).
#' @param add_tss_annotation Logical; whether to include TSS annotation columns in the output (default: FALSE).
#' @param verbose Logical; whether to print messages during execution (default: TRUE).
#'
#' @return A tibble or list containing RETSI, GETSI, and/or linked results depending on input and parameters.
#'
#' @export
run_spicey <- function(atac = NULL,
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
    stop("At least one of 'atac' or 'rna' data must be provided.")
  }

  if (!is.null(rna) && is.null(gene_id)) {
    stop("Parameter 'gene_id' must be provided when RNA data is supplied.")
  }

  # Compute RETSI
  if (!is.null(atac)) {
    retsi <- spicey_retsi(atac)
  }

  # Compute GETSI
  if (!is.null(rna)) {
    getsi <- spicey_getsi(rna, gene_id)
  }

  # If no annotation method, return based on inputs
  if (is.null(annot_method)) {
    if (!is.null(atac) && is.null(rna)) {
      message("🌶 SPICEY pipeline successfully completed")
      return(retsi)
    }
    if (is.null(atac) && !is.null(rna)) {
      message("🌶 SPICEY pipeline successfully completed")
      return(getsi)
    }
    if (!is.null(atac) && !is.null(rna)) {
      message("🌶 SPICEY pipeline successfully completed")
      return(list(retsi = retsi, getsi = getsi))
    }
  }

  # Annotation requires ATAC data
  if (is.null(atac)) {
    stop("ATAC data must be provided for region-to-gene annotation.")
  }

  if (!(annot_method %in% c("nearest", "coaccessibility"))) {
    stop("annot_method must be either 'nearest' or 'coaccessibility'.")
  }

  if (annot_method == "nearest") {
    retsi_annotated <- annotate_with_nearest(
      retsi = retsi,
      txdb = txdb,
      keep_mito = keep_mito,
      annot_dbi = annot_dbi,
      protein_coding_only = filter_protein_coding,
      verbose = verbose,
      add_tss_annotation = add_tss_annotation
    )
  }

  if (annot_method == "coaccessibility") {
    if (is.null(links)) {
      stop("Links must be provided when using coaccessibility annotation method.")
    }
    retsi_annotated <- annotate_with_coaccessibility(
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
  }

  # Now handle linking or returning annotated + getsi

  if (link_spicey_measures) {
    if (is.null(rna)) {
      stop("RNA data must be provided for linking RETSI to gene expression (GETSI).")
    }
    if (verbose) message("→ Computing GETSI (if not already computed)...")
    if (!exists("getsi")) {
      getsi <- spicey_getsi(rna, gene_id)
    }

    if (verbose) message("→ Linking RETSI and GETSI via ", annot_method, "...")

    combined <- NULL
    if (annot_method == "nearest") {
      combined <- link_spicey_nearest(
        retsi_annotated_nearest = retsi_annotated,
        getsi = getsi
      )
    } else {
      combined <- link_spicey_coaccessible(
        retsi_annotated_coacc = retsi_annotated,
        getsi = getsi
      )
    }
    message("🌶 SPICEY pipeline successfully completed")
    return(combined)
  } else {
    if (!is.null(rna)) {
      message("🌶 SPICEY pipeline successfully completed")
      return(list(retsi_annotated = retsi_annotated, getsi = getsi))
    } else {
      message("🌶 SPICEY pipeline successfully completed")
      return(retsi_annotated)
    }
  }
}
