#' Example single-cell ATAC-seq differential accessibility data
#'
#' A toy or example dataset representing single-cell ATAC-seq differential analysis results.
#' Each row corresponds to a chromatin peak assessed for differential accessibility across one or more cell types.
#' @docType data
#' @name atac
#' @usage data(atac)
#' @format A data frame or GRanges-like object with the following required columns:
#' \describe{
#'   \item{seqnames}{Chromosome name (e.g., "chr1")}
#'   \item{start}{Start coordinate of the peak}
#'   \item{end}{End coordinate of the peak}
#'   \item{avg_log2FC}{Average log2 fold-change in accessibility for the peak in the specific cell type}
#'   \item{p_val}{Raw or adjusted p-value from differential accessibility analysis}
#'   \item{cell_type}{Cell type label associated with the differential measurement}
#' }
#' Multiple rows may exist per peak, each for a different cell type.
#'
#' @source Simulated differential ATAC-seq data for SPICEY
"atac"


#' Example single cell RNA-seq differential expression data
#'
#' A toy or example dataset representing single-cell RNA-seq differential expression results.
#' Each row corresponds to a gene tested across one or more cell types.
#' @docType data
#' @name rna
#' @usage data(rna)
#' @format A data frame or GRanges-like object with the following required columns:
#' \describe{
#'   \item{gene_id}{Gene identifier (e.g., official gene symbol or Ensembl ID)}
#'   \item{avg_log2FC}{Average log2 fold-change in expression for the gene in the specific cell type}
#'   \item{p_val}{Raw or adjusted p-value from differential expression analysis}
#'   \item{cell_type}{Cell type label associated with the differential measurement}
#' }
#' Multiple rows may exist per gene, each for a different cell type.
#'
#' @source Simulated RNA-seq data for SPICEY
"rna"


#' Example Cicero co-accessibility links
#'
#' A toy or example dataset representing co-accessibility links inferred from tools such as Cicero.
#' These links are used in SPICEY for associating regulatory elements to target genes.
#' @docType data
#' @name cicero_links
#' @usage data(cicero_links)
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Peak1}{Genomic coordinate or identifier of the first peak}
#'   \item{Peak2}{Genomic coordinate or identifier of the second peak}
#'   \item{coaccess}{Co-accessibility score or correlation value between the two peaks}
#' }
#' Peaks must correspond exactly to those used in the ATAC-seq input.
#'
#' @source Simulated Cicero results on example ATAC-seq data
"cicero_links"
