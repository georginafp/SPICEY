
#' Link RETSI peaks with GETSI values from coaccessibility
#'
#' This function combines GETSI (Gene Expression Tissue Specificity Index)
#' data with RETSI (Regulatory Element Tissue Specificity Index) results,
#' matching on gene symbols and cell types. It renames entropy columns,
#' performs joins to merge relevant annotations, and resolves coordinate conflicts.
#'
#' @param re A data.frame or tibble containing RETSI results with at least the columns:
#'   - `genes_HPAP`: gene symbols
#'   - `cell_type`
#'   - `norm_entropy`
#'   - `seqnames`, `start`, `end`, `annotation`, `distanceToTSS`, `region`, `RETSI`
#'
#' @param getsi A data.frame or tibble containing GETSI results with at least the columns:
#'   - `symbol`: gene symbol
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
#' combined_df <- link_spicey_coaccessible(re, getsi)
#' }
#'
#' @export
link_spicey_coaccessible <- function(re, getsi) {
  message("→ Link SPICEY measures from coaccessibility...")

  getsi2add <- as.data.frame(getsi) |>
    dplyr::select(symbol, GETSI, cell_type, norm_entropy) |>
    dplyr::rename(GETSI_entropy = norm_entropy)

  result <- re |>
    data.frame() |>
    dplyr::rename(RETSI_entropy = norm_entropy) |>
    dplyr::right_join(getsi2add, by = c("genes_HPAP" = "symbol", "cell_type")) |>
    dplyr::select(seqnames, start, end, everything()) |>
    dplyr::left_join(
      getsi |>
        data.frame() |>
        dplyr::select(symbol, seqnames, start, end) |>
        dplyr::distinct(),
      by = c("genes_HPAP" = "symbol")
    ) |>
    dplyr::mutate(
      seqnames = dplyr::coalesce(seqnames.x, seqnames.y),
      start = dplyr::coalesce(start.x, start.y),
      end = dplyr::coalesce(end.x, end.y)
    ) |>
    dplyr::select(-seqnames.x, -seqnames.y, -start.x, -start.y, -end.x, -end.y) |>
    dplyr::select(seqnames, start, end, dplyr::everything()) |>
    regioneR::toGRanges()

  return(result)
}



#' Link RETSI peaks with GETSI values from nearest genes
#'
#' @param retsi GRanges with RETSI scores and nearest genes.
#' @param getsi GRanges with GETSI scores by gene and cell type.
#'
#' @return Annotated GRanges with RETSI and matched GETSI scores.
#' @export
link_spicey_nearest <- function(retsi, getsi) {
  message("→ Link SPICEY measures from nearest gene...")

  retsi_df <- as.data.frame(retsi) |>
    dplyr::select(seqnames, start, end, everything()) |>
    dplyr::rename(
      genes_nearest = nearestGeneSymbol,
      RETSI_entropy = norm_entropy
    )

  getsi_df <- as.data.frame(getsi) |>
    dplyr::select(symbol, GETSI, cell_type, norm_entropy) |>
    dplyr::rename(
      genes_nearest = symbol,
      GETSI_entropy = norm_entropy
    )

  annotated_df <- dplyr::left_join(retsi_df, getsi_df,
                                   by = c("genes_nearest", "cell_type")) |>
    regioneR::toGRanges()
  return(annotated_df)
}
