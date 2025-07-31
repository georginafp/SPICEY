#' Link RETSI regions to GETSI scores using gene-based association methods
#'
#' This function connects regulatory regions scored with RETSI
#' (Regulatory Element Tissue Specificity Index) to genes scored with
#' GETSI (Gene Expression Tissue Specificity Index), based on either:
#' \itemize{
#'   \item \code{"nearest"} — links to the closest gene using the \code{nearestGeneSymbol} column.
#'   \item \code{"coaccessibility"} — links based on co-accessibility networks using the \code{gene_coacc} column.}
#' @param retsi_annotated A \code{GRanges} or \code{data.frame} of regulatory regions
#' with computed RETSI scores. Must include at least:
#'   \itemize{
#'     \item \code{RETSI} — Regulatory Element Tissue Specificity Index score.
#'     \item \code{RETSI_entropy} — Entropy-based specificity score for accessibility.
#'     \item \code{nearestGeneSymbol} or \code{gene_coacc} — Gene name or ID used for linking
#'       depending on \code{method}:
#'       \itemize{
#'         \item If \code{method = "nearest"}, uses \code{nearestGeneSymbol}.
#'         \item If \code{method = "coaccessibility"}, uses \code{gene_coacc}}}
#' @param getsi A \code{data.frame} containing GETSI results, with at least the following columns:
#'   \describe{
#'     \item{\code{[gene_id]}}{Gene identifier column (e.g., gene symbol or Ensembl ID).
#'     The identifier format must match the one used in the selected linking method.}
#'     \item{\code{cell_type}}{Cell type or cluster name (used for merging).}
#'     \item{\code{GETSI}}{Gene Expression Tissue Specificity Index score.}
#'     \item{\code{GETSI_entropy}}{Entropy-based specificity score for expression.}}
#' @param method Character string specifying the linking strategy:
#'   \itemize{
#'     \item \code{"nearest"} — link peaks to the nearest annotated gene.
#'     \item \code{"coaccessibility"} — link peaks using co-accessibility associations.}
#'   The gene identifiers in \code{getsi} must be consistent with those used in the selected method
#'   (e.g., both should use gene symbols or both Ensembl IDs) to avoid mismatched links.
#' @return A \code{data.frame} combining RETSI and GETSI information, including:
#' \itemize{
#'   \item Regulatory region coordinates (\code{seqnames}, \code{start}, \code{end}).
#'   \item Cell type.
#'   \item Linked gene (\code{gene_linked}).
#'   \item RETSI and RETSI_entropy.
#'   \item GETSI and GETSI_entropy.}
#' @examples
#' data(atac)
#' data(rna)
#' data(cicero_links)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' getsi <- compute_spicey_index(rna=rna, gene_id = "gene_id")
#' retsi <- compute_spicey_index(atac=atac)
#'
#' retsi_gene_nearest <- annotate_with_nearest(retsi = retsi,
#'                                             txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'                                             annot_dbi = org.Hs.eg.db::org.Hs.eg.db)
#' spicey_nearest <- link_spicey(retsi_annotated = retsi_gene_nearest,
#'                               getsi = getsi,
#'                               method = "nearest")
#'
#' retsi_gene_coacc <- annotate_with_coaccessibility(links = cicero_links,
#'                                                   retsi = retsi,
#'                                                   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'                                                   annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
#'                                                   coaccess_cutoff_override = 0.25)
#' spicey_coacc <- link_spicey(
#'   retsi_annotated = retsi_gene_coacc,
#'   getsi = getsi,
#'   method = "coaccessibility")
#' @export
#' @importFrom dplyr rename right_join select coalesce
# link_spicey <- function(retsi_annotated,
#                         getsi,
#                         method = NULL) {
#   if (is.null(method)) {
#     stop("You must provide a value for 'method': either 'coaccessibility' or 'nearest'.")
#   }
#   method <- match.arg(method, choices = c("coaccessibility", "nearest"))
#   message("Linking SPICEY measures using method: ", method)
#   df <- retsi_annotated |> data.frame()
#   if (method == "nearest") {
#     df <- df |> dplyr::rename(gene_linked = nearestGeneSymbol)
#   } else if (method == "coaccessibility") {
#     df <- df |> dplyr::rename(gene_linked = gene_coacc)
#   }
#   result <- df |>
#     dplyr::right_join(
#       getsi |>
#         data.frame() |>
#         dplyr::select(gene_id, GETSI, cell_type, GETSI_entropy),
#       by = c("gene_linked" = "gene_id", "cell_type")) |>
#     dplyr::select(seqnames, start, end, everything())
#   return(result)
# }

