#' Link RETSI peaks with GETSI values based on gene linkage method
#'
#' This generalized function links RETSI (Regulatory Element Tissue Specificity Index) regions
#' with GETSI (Gene Expression Tissue Specificity Index) scores based on either coaccessibility
#' or nearest gene annotation.
#'
#' @param retsi_annotated A \code{data.frame} object containing RETSI results
#' including one of: `gene_coacc` (for coaccessibility) or `nearestGeneSymbol` (for nearest).
#' @param getsi A \code{data.frame} or tibble containing GETSI results with columns:
#'   - `gene_id`
#'   - `cell_type`
#'   - `GETSI`
#'   - `GETSI_entropy`
#' @param method Character; either `"coaccessibility"` or `"nearest"` to indicate
#' the gene linkage method used to associate RETSI regions with GETSI scores.
#' If `"coaccessibility"`, links are based on co-accessible genes (`gene_coacc`);
#' if `"nearest"`, links are based on the nearest gene symbol (`nearestGeneSymbol`).
#' @return A \code{data.frame} combining RETSI and GETSI specificity scores,
#' linked by gene identifier and including entropy values for each.
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
#'   method = "coaccessibility"
#' )
#' @export
#' @importFrom dplyr rename right_join select coalesce
link_spicey <- function(retsi_annotated,
                        getsi,
                        method = NULL) {
  if (is.null(method)) {
    stop("You must provide a value for 'method': either 'coaccessibility' or 'nearest'.")
  }

  method <- match.arg(method, choices = c("coaccessibility", "nearest"))
  message("Linking SPICEY measures using method: ", method)

  df <- retsi_annotated |> data.frame()

  # Rename gene column according to method
  if (method == "nearest") {
    df <- df |> dplyr::rename(gene_linked = nearestGeneSymbol)
  } else if (method == "coaccessibility") {
    df <- df |> dplyr::rename(gene_linked = gene_coacc)
  }

  result <- df |>
    dplyr::right_join(
      getsi |>
        data.frame() |>
        dplyr::select(gene_id, GETSI, cell_type, GETSI_entropy),
      by = c("gene_linked" = "gene_id", "cell_type")
    ) |>
    dplyr::select(seqnames, start, end, everything())

  return(result)
}

