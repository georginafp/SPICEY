#' Example single-cell ATAC-seq differential accessibility data
#'
#' A simulated example dataset representing single-cell ATAC-seq
#' differential accessibility results.
#' Each row corresponds to a chromatin accessibility peak tested for
#' differential accessibility across one or more cell types.
#' The dataset is formatted as a list of \code{GRanges} objects or
#' data frames (convertible to \code{GRanges}),  where each element represents
#' differential accessibility statistics for a specific cell type.
#' Multiple rows may exist per peak, each representing results in a different cell type.
#' @docType data
#' @name atac
#' @usage data(atac)
#' @format A data frame or \code{GRanges}-like object with the following required columns:
#' \describe{
#'   \item{seqnames}{Chromosome name (e.g., "chr1")}
#'   \item{start}{Start coordinate of the peak}
#'   \item{end}{End coordinate of the peak}
#'   \item{avg_log2FC}{Average log2 fold-change of accessibility for the peak in the specific cell type}
#'   \item{p_val}{Raw or adjusted p-value from differential accessibility testing}
#'   \item{cell_type}{Cell type or cluster label associated with each measurement}
#' }
#'
#' @source Simulated single-cell ATAC-seq differential accessibility data used in SPICEY examples.
"atac"


#' Example single-cell RNA-seq differential expression data
#'
#' A simulated example dataset representing single-cell RNA-seq
#' differential expression results.
#' Each row corresponds to a gene tested across one or more cell types.
#' The dataset is formatted as a list of \code{GRanges} objects or
#' data frames (convertible to \code{GRanges}), where each element contains
#' differential expression statistics for a specific cell type.
#' Multiple rows may exist per gene, each representing results in a different cell type.
#' @docType data
#' @name rna
#' @usage data(rna)
#' @format A data frame or \code{GRanges}-like object with the following required columns:
#' \describe{
#'   \item{gene_id}{Gene identifier (e.g., official gene symbol or Ensembl ID)}
#'   \item{avg_log2FC}{Average log2 fold-change of expression for the gene in the specific cell type}
#'   \item{p_val}{Raw or adjusted p-value from differential expression testing}
#'   \item{cell_type}{Cell type or cluster label associated with each measurement}
#' }
#' @source Simulated single-cell RNA-seq differential expression data used in SPICEY examples.
"rna"


#' Example Cicero co-accessibility links
#'
#' A simulated example dataset of co-accessibility links inferred
#' from single-cell ATAC-seq data using tools such as Cicero or Signac's LinkPeaks().
#' These links support integrative analysis by associating regulatory elements
#' with putative target genes. The peaks referenced here must exactly match
#' those in the ATAC-seq differential accessibility dataset.
#' @docType data
#' @name cicero_links
#' @usage data(cicero_links)
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Peak1}{Genomic coordinate or peak identifier for the first peak in the pair}
#'   \item{Peak2}{Genomic coordinate or peak identifier for the second peak in the pair}
#'   \item{coaccess}{Co-accessibility score or correlation value quantifying the linkage}
#' }
#' @source Simulated Cicero co-accessibility results on example single-cell ATAC-seq data.
"cicero_links"
