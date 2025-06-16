#' Retrieve gene annotations from Ensembl via biomaRt
#'
#' Downloads gene annotation data from Ensembl using biomaRt, processes it,
#' and returns a list with both a data.frame and a filtered GRanges object
#' containing protein-coding genes. This is useful for RNA-seq and epigenomic
#' pipelines.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{df}{A \code{data.frame} with gene annotation information for RNA-seq analysis.}
#'   \item{gr}{A \code{GRanges} object with protein-coding genes for epigenomic analyses.}
#' }
#'
#' @import biomaRt
#' @import regioneR
#' @import GenomeInfoDb
#' @importFrom plyranges filter
#' @export
biomart_genes <- function() {
  mart <- biomaRt::useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    host = "https://www.ensembl.org",
    path = "/biomart/martservice",
    dataset = "hsapiens_gene_ensembl"
  )

  genes <- biomaRt::getBM(
    attributes = c(
      "chromosome_name",
      "start_position", "end_position", "strand",
      "ensembl_gene_id", "external_gene_name", "gene_biotype",
      "percentage_gene_gc_content"
    ),
    useCache = TRUE,
    mart = mart
  )

  message("Got it!")

  ## Convert strand from numeric to factor with "+" and "-"
  genes$strand <- ifelse(genes$strand == 1, "+",
                         ifelse(genes$strand == -1, "-", "*"))

  ## Convert to GRanges
  genes_gr <- regioneR::toGRanges(genes)
  strand(genes_gr) <- genes_gr$strand
  mcols(genes_gr) <- mcols(genes_gr)[, setdiff(colnames(mcols(genes_gr)), "strand"), drop=FALSE]

  ## Keep only standard chromosomes except MT (mitochondrial)
  genes_gr <- GenomeInfoDb::keepStandardChromosomes(genes_gr, pruning.mode = "coarse")
  genes_gr <- GenomeInfoDb::dropSeqlevels(genes_gr, "MT", pruning.mode = "coarse")
  GenomeInfoDb::seqlevelsStyle(genes_gr) <- "UCSC"

  ## Filter protein coding genes
  genes_gr <- plyranges::filter(genes_gr, gene_biotype == "protein_coding")

  ## Select columns for output GRanges
  genes_gr <- genes_gr[, c("ensembl_gene_id", "external_gene_name")]

  ## Return list
  list(df = genes, gr = genes_gr)
}


#' Read input data from a file path
#'
#' Reads data from the specified file path. Supports RDS, CSV, TSV, and TXT formats.
#' Returns a data.frame or a GRanges object if RDS contains GRanges.
#'
#' @param path Character string of the file path to read.
#'
#' @return The contents of the file, either a \code{data.frame}, \code{GRanges}, or other R object saved in the RDS.
#'
#' @importFrom tools file_ext
#' @importFrom data.table fread
#' @export
read_input_file <- function(path) {
  ext <- tolower(file_ext(path))
  if (ext == "rds") {
    message("Reading RDS: ", path)
    readRDS(path)
  } else if (ext %in% c("csv")) {
    message("Reading CSV: ", path)
    fread(path, data.table = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    message("Reading TSV: ", path)
    fread(path, sep = "\t", data.table = FALSE)
  } else {
    stop("Unsupported file type: ", ext)
  }
}
