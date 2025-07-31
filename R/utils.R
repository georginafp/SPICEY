utils::globalVariables(c(
  "peak", "seqnames", "start", "end", "distanceToTSS", "annotation",
  "avg_FC", "max_FC", "norm_FC", "weight", "entropy", "score",
  "norm_entropy", ".data", "prob", "entropy_component", "name", "ENTREZID",
  "GENETYPE", "nearestGeneSymbol", "gene_coacc", "gene_id", "GETSI",
  "cell_type", "GETSI_entropy", "avg_log2FC", "p_val", "p_val_adj", "everything"
))


#' Get gene promoters with optional filtering for protein-coding genes
#'
#' @param txdb A TxDb object.
#' @param annot_dbi AnnotationDbi object (e.g., org.Hs.eg.db).
#' @param upstream Upstream window from TSS.
#' @param downstream Downstream window from TSS.
#' @param protein_coding_only Logical, whether to filter to protein-coding genes.
#' @return GRanges of promoters annotated with gene symbol.
get_promoters <- function(txdb,
                          annot_dbi,
                          upstream = 2000,
                          downstream = 2000,
                          protein_coding_only = TRUE) {

  proms <- GenomicFeatures::promoters(
    GenomicFeatures::genes(txdb),
    upstream = upstream,
    downstream = downstream)

  symbols <- AnnotationDbi::mapIds(
    annot_dbi,
    keys = names(proms),
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first")

  if (protein_coding_only) {
    gene_types <- AnnotationDbi::mapIds(
      annot_dbi,
      keys = names(proms),
      column = "GENETYPE",
      keytype = "ENTREZID",
      multiVals = "first")
    keep <- which(gene_types == "protein-coding")
    proms <- proms[keep]
    symbols <- symbols[keep]
  }

  proms$gene_id <- symbols
  return(proms)
}



#' Convert GRanges to String Representation
#'
#' Converts a \code{GRanges} object into a character vector where each element
#' represents a genomic range in the format "chr-start-end".
#'
#' @param gr A \code{GRanges} object.
#'
#' @return A character vector of genomic ranges in "chr-start-end" format.
#'
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = c("chr1", "chr2"),
#'               ranges = IRanges(start = c(100, 200), end = c(150, 250)))
#' granges_to_string(gr)
#' # [1] "chr1-100-150" "chr2-200-250"
#'
#' @export
granges_to_string <- function(gr) {
  if (!inherits(gr, "GRanges")) {
    stop("Input must be a GRanges object")
  }
  string <- paste0(as.character(GenomicRanges::seqnames(gr)),
                   "-", GenomicRanges::start(gr),
                   "-", GenomicRanges::end(gr))
  return(string)
}




#' Parse inputs from the different accepted types into data.frames
#' @importFrom GenomicRanges mcols
.parse_input_diff <- function(input) {
  if(is(input, "list")) {
    if(is.null(names(input)))  stop("If your differential regions are in a list it should be named with cell types")

    if(is(input[[1]], "GRanges")) {
      input_df <- lapply(input, function(x) data.frame(mcols(x)))

      input <- dplyr::bind_rows(input,
                                .id = "cell_type")
    } else if (is(input[[1]], "data.frame")) {
      input <- dplyr::bind_rows(input,
                                .id = "cell_type")
    }
  } else if (is(input, "GRangesList")) {
    input <- dplyr::bind_rows(lapply(input, function(x) data.frame(mcols(x))),
                              .id = "cell_type")
  }

  return(input)
}
