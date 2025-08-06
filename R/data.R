#' Example single-cell ATAC-seq differential accessibility data
#'
#' A toy example dataset representing single-cell ATAC-seq differential accessibility results.
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
#'   \item{region_id}{Unique identifier of the region (e.g., \code{chr1-5000-5800}).}
#'   \item{avg_log2FC}{Average log2 fold-change of accessibility for the peak in the specific cell type}
#'   \item{p_val_adj}{Adjusted p-value (e.g., \code{FDR-corrected})}
#'   \item{cell_type}{Cell type or cluster label associated with each measurement (e.g., \code{Acinar})}
#' }
#' @source
#' Precomputed using \code{FindMarkers()} (Wilcoxon test, via Presto if available) on control samples
#' from the Human Pancreas Analysis Program (HPAP), using paired snATAC-seq and snRNA-seq data
#' from three non-diabetic human donors.
"atac"


#' Example single-cell RNA-seq differential expression data
#'
#' A toy example dataset representing single-cell RNA-seq
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
#'   \item{gene_id}{Identifier of the gene. Must be an official gene symbol (e.g., \code{GAPDH}).
#'   \item{avg_log2FC}{Average log2 fold-change of expression for the gene in the specific cell type}
#'   \item{p_val_adj}{Adjusted p-value (e.g., \code{FDR-corrected})}
#'    \item{cell_type}{Cell type or cluster label associated with each measurement (e.g., \code{Acinar})}
#' }
#' @source
#' Precomputed using \code{FindMarkers()} (Wilcoxon test, via Presto if available) on control samples
#' from the Human Pancreas Analysis Program (HPAP), using paired snATAC-seq and scRNA-seq data
#' from three non-diabetic human donors.
"rna"


#' Example Cicero co-accessibility links
#'
#' A toy example dataset of co-accessibility links inferred
#' from single-cell ATAC-seq data using tools such as Cicero or Signac's LinkPeaks().
#' These links support integrative analysis by associating regulatory elements
#' with putative target genes. The peaks referenced here must exactly match
#' those in the ATAC-seq differential accessibility dataset.
#' @docType data
#' @name cicero_links
#' @usage data(cicero_links)
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Peak1}{Genomic coordinate or peak identifier for the first peak in the pair (e.g., \code{chr1-110209621-110211746})}
#'   \item{Peak2}{Genomic coordinate or peak identifier for the second peak in the pair (e.g., \code{chr1-110209621-110211746})}
#'   \item{coaccess}{Co-accessibility score or correlation value quantifying the linkage}
#' }
#' @source
#' Cicero co-accessibility links were computed from UMAP-reduced snATAC-seq data
#' (HPAP, control donors) using \code{run_cicero()} with chromosome sizes from hg38.
#' Input data matched the peaks in the provided ATAC dataset.
"cicero_links"
