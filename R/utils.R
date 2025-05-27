##..............................................................................
# BiomaRt gene data                                                         ####
##..............................................................................
#'
#' @return List with two elements:
#'         - `df`. Data.frame with gene information to use for the RNA-seq
#'         pipeline.
#'         - `gr`. GRanges with protein coding genes to use for the epigenome
#'         pipeline and analyses.
biomart_genes <- function() {
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           host = "https://www.ensembl.org",
                           path = "/biomart/martservice",
                           dataset = "hsapiens_gene_ensembl")
  genes <- biomaRt::getBM(attributes = c("chromosome_name",
                                         "start_position", "end_position", "strand",
                                         "ensembl_gene_id", "external_gene_name", "gene_biotype",
                                         "percentage_gene_gc_content"),
                          useCache = TRUE, mart = mart)
  message("Got it!")

  ## Make GRanges
  genes$strand[genes$strand == -1] <- "-"
  genes$strand[genes$strand == 1] <- "+"

  genes_gr <- regioneR::toGRanges(genes)
  strand(genes_gr) <- genes_gr$strand
  mcols(genes_gr) <- mcols(genes_gr)[, -1]
  genes_gr <- GenomeInfoDb::keepStandardChromosomes(genes_gr, pruning.mode = "coarse")
  genes_gr <- GenomeInfoDb::dropSeqlevels(genes_gr, "MT", pruning.mode = "coarse")
  GenomeInfoDb::seqlevelsStyle(genes_gr) <- "UCSC"

  genes_gr <- genes_gr %>%
    plyranges::filter(gene_biotype == "protein_coding")

  genes_gr <- genes_gr[,c("ensembl_gene_id", "external_gene_name")]

  ## Make list with both
  genes_list <- list(df=genes,
                     gr=genes_gr)
  return(genes_list)

}


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
