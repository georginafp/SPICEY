#' Run full SPICEY pipeline: compute RETSI and GETSI, annotate regulatory elements, and link measures
#'
#' This function performs the complete SPICEY pipeline workflow:
#' 1. Computes Regulatory Element Tissue Specificity Index (RETSI) from scATAC-seq DE results.
#' 2. Computes Gene Expression Tissue Specificity Index (GETSI) from scRNA-seq DE results.
#' 3. Annotates regulatory elements (RE) with nearest genes.
#' 4. Optionally annotates RE with co-accessible target genes if `link_method = "coaccessible"`.
#' 5. Links RETSI and GETSI scores either by nearest gene or co-accessibility links.
#'
#' @param atac_data A data.frame or tibble containing differential accessibility results used to compute RETSI.
#' @param rna_data A data.frame or tibble containing differential expression results used to compute GETSI.
#' @param re_gr A GRanges object representing regulatory elements to be annotated and linked.
#' @param gene_gr A GRanges object representing genes used for nearest or co-accessibility annotation.
#' @param link_method Character string specifying linking method. One of `"nearest"` (default) or `"coaccessible"`.
#' @param coaccess_links Optional. A data.frame or GRangesList representing co-accessibility links. Required if `link_method = "coaccessible"`.
#'
#' @return A list containing:
#' \describe{
#'   \item{retsi}{Data frame with computed RETSI scores.}
#'   \item{getsi}{Data frame with computed GETSI scores.}
#'   \item{re_gr_annotated}{GRanges object with regulatory elements annotated with nearest and/or co-accessible genes.}
#'   \item{linked}{Data frame linking RETSI and GETSI scores according to chosen linking method.}
#' }
#'
#' @examples
#' \dontrun{
#' # Run with nearest gene linking (default)
#' results <- run_spicey_full_pipeline(atac_data, rna_data, re_gr, gene_gr)
#'
#' # Run with co-accessibility linking
#' results <- run_spicey_full_pipeline(
#'   atac_data, rna_data, re_gr, gene_gr,
#'   link_method = "coaccessible",
#'   coaccess_links = coaccessibility_links
#' )
#' }
#'
#' @export
run_spicey <- function(
    atac_data,
    rna_data,
    re_gr,
    gene_gr,
    link_method = c("nearest", "coaccessible"),
    coaccess_links = NULL
) {
  link_method <- match.arg(link_method)

  # 1. Compute RETSI
  retsi_df <- spicey_retsi(atac_data)

  # 2. Compute GETSI
  getsi_df <- spicey_getsi(rna_data)

  # 3. Annotate regulatory elements with nearest gene
  re_gr <- annotate_with_nearest(re_gr, gene_gr)

  # 4. Annotate with co-accessible target genes if requested
  if (link_method == "coaccessible") {
    if (is.null(coaccess_links)) {
      stop("coaccess_links must be provided if link_method = 'coaccessible'")
    }
    re_gr <- annotate_with_coaccessibility(re_gr, gene_gr, coaccess_links)
  }

  # 5 & 6. Link spicey measures based on chosen method
  if (link_method == "nearest") {
    linked <- link_spicey_nearest(retsi_df, getsi_df, re_gr)
  } else {
    linked <- link_spicey_coaccessibility(retsi_df, getsi_df, re_gr)
  }

  # Return results list
  list(
    retsi = retsi_df,
    getsi = getsi_df,
    re_gr_annotated = re_gr,
    linked = linked
  )
}
