#' Run full SPICEY pipeline: compute RETSI and GETSI, annotate regulatory elements, and link measures
#'
#' This function performs the complete SPICEY pipeline workflow:
#' 1. Computes Regulatory Element Tissue Specificity Index (RETSI) from scATAC-seq differential accessibility results.
#' 2. Computes Gene Expression Tissue Specificity Index (GETSI) from scRNA-seq differential expression results.
#' 3. Annotates regulatory elements (RE) with nearest genes.
#' 4. Optionally annotates RE with co-accessible target genes if `link_method = "coaccessibility"`.
#' 5. Links RETSI and GETSI scores either by nearest gene or co-accessibility links.
#'
#' @param atac A \code{GRanges} or \code{data.frame} convertible to \code{GRanges} containing differential accessibility results. Must include \code{seqnames}, \code{start}, \code{end}, \code{avg_log2FC}, \code{p_val}, and \code{cell_type} columns.
#' @param rna A \code{GRanges} or \code{data.frame} convertible to \code{GRanges} containing differential expression results. Must include \code{seqnames}, \code{start}, \code{end}, \code{avg_log2FC}, \code{p_val}, and \code{cell_type} columns.
#' @param link_method Character string specifying the linking method. One of \code{"nearest"} or \code{"coaccessibility"}. If \code{NULL} (default), linking is not performed.
#' @param links Optional. A \code{data.frame} or \code{GRangesList} of co-accessibility links. Required if \code{link_method = "coaccessibility"}.
#' @param coaccess_cutoff_override Numeric scalar. Co-accessibility score cutoff used to filter links. Default is \code{0.25}.
#' @param filter_promoter_distal Logical. Whether to filter links to keep only promoter-distal interactions when using co-accessibility. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{retsi}{Data frame with computed RETSI scores if \code{atac} provided.}
#'   \item{getsi}{Data frame with computed GETSI scores if \code{rna} provided.}
#'   \item{linking}{Data frame linking RETSI and GETSI scores if \code{link_method} specified.}
#' }
#'
#' @examples
#' \dontrun{
#' # Compute GETSI only
#' results <- run_spicey_pipeline(rna = rna_data)
#'
#' # Compute RETSI only
#' results <- run_spicey_pipeline(atac = atac_data)
#'
#' # Compute both GETSI and RETSI without linking
#' results <- run_spicey_pipeline(atac = atac_data, rna = rna_data)
#'
#' # Compute both and link with nearest gene method
#' results <- run_spicey_pipeline(atac = atac_data, rna = rna_data, link_method = "nearest")
#'
#' # Compute both and link with co-accessibility method
#' results <- run_spicey_pipeline(
#'   atac = atac_data,
#'   rna = rna_data,
#'   link_method = "coaccessibility",
#'   links = coaccessibility_links,
#'   coaccess_cutoff_override = 0.3,
#'   filter_promoter_distal = TRUE
#' )
#' }
#'
#' @export
run_spicey <- function(atac = NULL,
                                rna = NULL,
                                link_method = NULL,
                                links = NULL,
                                coaccess_cutoff_override = 0.25,
                                filter_promoter_distal = TRUE) {
  # Load TxDb for annotation (adjust as needed)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

  # Check inputs
  if (is.null(atac) && is.null(rna)) {
    stop("At least one of 'atac' or 'rna' data must be provided.")
  }

  # Initialize output list
  results <- list()

  # Compute RETSI if atac provided
  if (!is.null(atac)) {
    message("Computing RETSI...")
    retsi <- spicey_retsi(atac)
    results$retsi <- retsi
  }

  # Compute GETSI if rna provided
  if (!is.null(rna)) {
    message("Computing GETSI...")
    getsi <- spicey_getsi(rna)
    results$getsi <- getsi
  }

  # If both atac and rna provided, and no link_method, just return both indices
  if (!is.null(atac) && !is.null(rna) && is.null(link_method)) {
    return(results)
  }

  # If linking requested, check required inputs and perform linking
  if (!is.null(link_method)) {
    if (is.null(atac) || is.null(rna)) {
      stop("Both 'atac' and 'rna' data must be provided for linking.")
    }
    if (!(link_method %in% c("nearest", "coaccessibility"))) {
      stop("link_method must be either 'nearest' or 'coaccessibility'")
    }

    if (link_method == "nearest") {
      message("Linking regulatory elements to genes using nearest gene method...")
      # Annotate RE with nearest genes
      retsi_gene_nearest <- annotate_with_nearest(results$retsi)
      # Link spicey using nearest method
      spicey_nearest <- link_spicey_nearest(retsi_gene_nearest, results$getsi)
      results$linking <- spicey_nearest
    }

    if (link_method == "coaccessibility") {
      if (is.null(links)) {
        stop("Argument 'links' must be provided when link_method = 'coaccessibility'")
      }
      message("Filtering and annotating coaccessibility links...")

      # Annotate links with CCANs and filter promoter-distal if requested
      annotated_links <- annotate_links_with_ccans(
        links = links,
        coaccess_cutoff_override = coaccess_cutoff_override,
        filter_promoter_distal = filter_promoter_distal,
        txdb = txdb
      )

      message("Linking regulatory elements to genes using coaccessibility method...")
      spicey_coacc <- link_spicey_coaccessible(
        re = results$retsi,
        getsi = results$getsi,
        links = annotated_links,
        txdb = txdb,
        name_links = "coacc"
      )
      results$linking <- spicey_coacc
    }
  }

  return(results)
}

