utils::globalVariables(c(
  "region_id", "seqnames", "start", "end", "distanceToTSS", "annotation",
  "avg_FC", "max_FC", "norm_FC", "weight", "entropy", "score",
  "norm_entropy", ".data", "prob", "entropy_component", "name", "ENTREZID",
  "GENETYPE", "gene_id", "gene_coacc", "gene_id", "GETSI",
  "Peak1", "Peak2", "coaccess", "gene_name1", "gene_name2",
  "peak", "queryHits", "subjectHits", "seqlevels", "is", "makeGRangesFromDataFrame",
  "cell_type", "GETSI_entropy", "avg_log2FC", "p_val", "p_val_adj", "GETSI_z",
  "RETSI", "RETSI_z", "combined_score", "max_cell_type", "max_score", "z_score", "everything"
))

#' Parses input data of various types (e.g., named lists of \code{GRanges}
#' or \code{data.frame}, or a \code{GRangesList}) into a single tidy
#' \code{data.frame}, with a \code{cell_type} column.
#' @param input An object representing differential results, such as:
#'   \itemize{
#'     \item A named list of \code{GRanges} objects.
#'     \item A named list of \code{data.frame}s.
#'     \item A \code{GRangesList}.
#'   }
#' @return A \code{data.frame} combining all elements, with an added \code{cell_type} column indicating the source.
#' @importFrom GenomicRanges mcols
.parse_input_diff <- function(input) {
  if (is(input, "list")) {
    if (is.null(names(input))) stop("If your differential regions are in a list it should be named with cell types")
    if (is(input[[1]], "GRanges")) {
      input_df <- lapply(input, function(x) {
        df <- data.frame(mcols(x))
        if ("cell_type" %in% colnames(df)) {
          df$cell_type <- NULL
        }
        df
      })
      input <- dplyr::bind_rows(input_df, .id = "cell_type")
    } else if (is(input[[1]], "data.frame")) {
      input <- dplyr::bind_rows(input, .id = "cell_type")
    }
  } else if (is(input, "GRangesList")) {
    input <- dplyr::bind_rows(lapply(input, function(x) {
      df <- data.frame(mcols(x))
      if ("cell_type" %in% colnames(df)) {
        df$cell_type <- NULL
      }
      df
    }), .id = "cell_type")
  }
  return(input)
}
