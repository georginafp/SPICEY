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




#' Extract promoter regions of protein-coding genes
#'
#' This function extracts ±2kb promoter regions around transcription start sites
#' from a `TxDb` object, restricts to standard chromosomes, and annotates each
#' region with `GENEID` and gene `symbol`. Only protein-coding genes are retained.
#'
#' @param txdb A `TxDb` object (e.g., `TxDb.Hsapiens.UCSC.hg38.knownGene`)
#'
#' @return A `GRanges` object with promoter regions annotated with `GENEID` and `symbol`
#' @export
get_promoters_protein_coding <- function(txdb) {
  # Extract transcripts and define ±2kb promoters
  proms <- GenomicFeatures::transcripts(txdb, columns = c("GENEID")) |>
    GenomicRanges::promoters(upstream = 2000, downstream = 2000) |>
    GenomeInfoDb::keepSeqlevels(paste0("chr", c(1:22, "X", "Y")),
                                pruning.mode = "coarse")

  # Replace empty GENEID with NA
  proms$GENEID <- sapply(proms$GENEID, function(x) ifelse(length(x) == 0, NA, x))

  # Add gene symbols from org.Hs.eg.db
  proms$symbol <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                        keys = proms$GENEID,
                                        columns = "SYMBOL",
                                        keytype = "ENTREZID",
                                        multiVals = "first")$SYMBOL

  # Identify protein-coding genes
  pc <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                              keys = proms$GENEID,
                              columns = "GENETYPE",
                              keytype = "ENTREZID",
                              multiVals = "first") |>
    dplyr::filter(!is.na(ENTREZID), !is.na(GENETYPE), GENETYPE == "protein-coding") |>
    dplyr::pull(ENTREZID)

  # Filter for protein-coding gene promoters
  proms <- proms[which(proms$GENEID %in% pc)]

  return(proms)
}
