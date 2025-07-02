
#' Link RETSI peaks with GETSI values from coaccessibility
#'
#' This function combines GETSI (Gene Expression Tissue Specificity Index)
#' data with RETSI (Regulatory Element Tissue Specificity Index) results,
#' matching on gene symbols and cell types. It renames entropy columns,
#' performs joins to merge relevant annotations, and resolves coordinate conflicts.
#'
#' @param re A data.frame or tibble containing RETSI results with at least the columns:
#'   - `gene_coacc`: gene symbols of the coaccessible gene to that region
#'   - `cell_type`
#'   - `norm_entropy`
#'   - `seqnames`, `start`, `end`, `annotation`, `distanceToTSS`, `region`, `RETSI`
#'
#' @param getsi A data.frame or tibble containing GETSI results with at least the columns:
#'   - `gene_id`: gene gene_id
#'   - `cell_type`
#'   - `GETSI`
#'   - `norm_entropy`
#'   - `seqnames`, `start`, `end`
#'
#' @return A data.frame with combined RETSI and GETSI information including entropy values,
#' gene annotations, and harmonized genomic coordinates.
#'
#' @import dplyr
#' @importFrom dplyr coalesce
#'
#' @examples
#' \dontrun{
#' combined_df <- link_spicey_coaccessible(retsi_annotated_coacc, getsi)
#' }
#'
#' @export
link_spicey_coaccessible <- function(retsi_annotated_coacc, getsi) {
  message("→ Link SPICEY measures from coaccessibility...")
  result <- retsi_annotated_coacc |>
    data.frame() |>
    dplyr::rename(RETSI_entropy = norm_entropy) |>
    dplyr::right_join(getsi |>
                        data.frame() |>
                        dplyr::select(gene_id, GETSI, cell_type, norm_entropy) |>
                        dplyr::rename(GETSI_entropy = norm_entropy),
                      by = c(gene_coacc = "gene_id", "cell_type")) |>
    dplyr::select(seqnames, start, end, everything()) |>
    data.frame()

  return(result)
}



#' Link RETSI peaks with GETSI values from nearest genes
#'
#' @param retsi_annotated_nearest GRanges with RETSI scores and nearest genes.
#' @param getsi GRanges with GETSI scores by gene and cell type.
#'
#' @return Annotated GRanges with RETSI and matched GETSI scores.
#' @export
link_spicey_nearest <- function(retsi_annotated_nearest, getsi) {
  message("→ Link SPICEY measures from nearest gene...")

  result <- retsi_annotated_nearest |>
    data.frame() |>
    dplyr::rename(gene_nearest = nearestGeneSymbol,
                  RETSI_entropy = norm_entropy) |>
    dplyr::right_join(getsi |>
                        data.frame() |>
                        dplyr::select(gene_id, GETSI, cell_type, norm_entropy) |>
                        dplyr::rename(GETSI_entropy = norm_entropy),
                      by = c(gene_nearest = "gene_id", "cell_type")) |>
    dplyr::select(seqnames, start, end, everything()) |>
    data.frame()

    return(result)
}
