#' SPICEY for tissue specificity analysis
#'
#' Computes tissue-specificity scores from differential accessibility (RETSI)
#' and/or gene expression (GETSI) data obtained from single cell experiments. Supports:
#' \itemize{
#'   \item RETSI calculation from differential accessibility data in different cell types/clusters (scATAC-seq).
#'   \item GETSI calculation from differential expression data in different cell types/clusters (scRNA-seq).
#'   \item Optional integration of RETSI and GETSI scores by linking gene associations (see \code{\link{link_atac_to_genes}}).
#' }
#' @param rna Either a single \code{data.frame} or a named list of 
#'   \code{data.frame}s or \code{GRanges} where each element corresponds 
#'   to a cell type. It should contain differential expression results, 
#'   with required columns:
#'   \describe{
#'     \item{gene_id}{Identifier of the gene (e.g., gene symbol, Ensembl ID). 
#'          The name of this column should be provided in argument \code{gene_id}}
#'     \item{avg_log2FC}{Average log2 fold-change for the gene in that cell type.}
#'     \item{p_val}{Raw p-value for the differential test.}
#'     \item{p_val_adj}{Adjusted p-value (e.g., FDR-corrected).}
#'     \item{cell_type}{Cell type or cluster label. Only necessary when input is 
#'           a single \code{data.frame}. If input is a list, it will be generated 
#'           from list names}.
#'   }
#'   Note that the same gene may appear multiple times across cell types.
#' @param atac Either a single \code{data.frame} or a named list of 
#'   \code{data.frame}s or \code{GRanges} where each element corresponds 
#'   to a cell type.  It should contain differential chromatin accessibility 
#'   results with required columns:
#'   \describe{
#'     \item{region_id}{Unique identifier of the region. The name of this column 
#'           should be provided in argument \code{region_id}.}
#'     \item{avg_log2FC}{Average log2 fold-change for accessibility in that cell type.}
#'     \item{p_val}{Raw p-value for the differential test.}
#'     \item{p_val_adj}{Adjusted p-value (e.g., FDR).}
#'     \item{cell_type}{Cell type or cluster label. Only necessary when input is 
#'           a single \code{data.frame}. If input is a list, it will be generated 
#'           from list names}.
#'   Note that the same region may appear multiple times across cell types.
#'   }
#' @param gene_id A character string specifying the column name in each list element
#'   that contains the gene identifiers.
#' @param region_id A character string specifying the column name in each list element
#'   that contains the accessible region identifiers.
#' @param annotation A data.frame linking \code{gene_id} to \code{region_id}. They should have 
#'   the same names provided in the respective parameters. This can be provided 
#'   by the user or generated using the function \code{\link{link_atac_to_genes}}.
#' @param verbose Logical; print messages (default TRUE).
#' @return Depending on inputs, returns RETSI and/or GETSI data frames, optionally linked and annotated.
#' @examples
#' data(atac)
#' retsi <- SPICEY(atac=atac, region_id="region_id")
#' head(retsi)
#'
#' data(rna)
#' getsi <- SPICEY(rna=rna, gene_id="gene_id")
#' head(getsi)
#' 
#' both <-  SPICEY(rna = rna, gene_id = "gene_id",
#'                 atac=atac, region_id="region_id")
#' lapply(both, head)
#' TODO: Include example with annotation
#' @export
SPICEY <- function(atac = NULL,
                   rna = NULL,
                   gene_id = NULL,
                   region_id = NULL,
                   annotation = NULL,
                   verbose = TRUE) {

  if (is.null(atac) && is.null(rna)) {
    stop("Provide at least one of 'atac' or 'rna'.")
  }

  if (!is.null(rna) && is.null(gene_id)) {
    stop("'gene_id' is required when RNA data is supplied.")
  }
  
  if (!is.null(atac) && is.null(region_id)) {
    stop("'region_id' is required when ATAC data is supplied.")
  }
  
  #CHANGE: Make SPICEY deal with whether the input is atac or rna. Input one at a 
  # time to compute_spicey_index --> For that function it does not matter whether
  # the data is atac or RNA
  
  # 1) Compute GETSI if RNA available
  if (!is.null(rna)) {
    message("Computing GETSI & entropy...")
    rna <- .parse_input_diff(rna)
    getsi <- compute_spicey_index(diff = rna, id = gene_id) |>
      dplyr::rename(GETSI = score,
                    GETSI_entropy = norm_entropy)
  } else {
    getsi <- NULL
  }
  
  # 2) Compute RETSI if ATAC available
  if (!is.null(atac)) {
    message("Computing RETSI & entropy...")
    atac <- .parse_input_diff(atac)
    retsi <- compute_spicey_index(atac, id = region_id) |>
      dplyr::rename(RETSI = score,
                    RETSI_entropy = norm_entropy)
  } else {
    retsi <- NULL
  }
  
  results <- list(RETSI = retsi, GETSI = getsi)
  results <- results[!sapply(results, is.null)]

  # 3) If no annotation, return these results directly (list if both are 
  # calculated and data.frame if only one of them is). If not, link measures
  if(is.null(annotation) & length(results) > 1) {
    # Return list if both RETSI and GETSI calculated (but no annotation)
    message("SPICEY pipeline successfully completed")
    return(results)
  } else if(is.null(annotation) & length(results) == 1) {
    # Return data.frame if only one of GETSI or RETSI were calculates
    message("SPICEY pipeline successfully completed")
    return(results[[1]])
  } else {
    # If annotation is provided, linke RETSI and GETSI
    message("Linking RETSI and GETSI using provided annotation...")
    combined <- link_spicey(retsi = retsi,
                            region_id = region_id,
                            getsi = getsi,
                            gene_id = gene_id,
                            annotation = annotation)
    # Return RE-gene links with RETSI and GETSI
    results$linked <- combined
    
    message("SPICEY pipeline successfully completed")
    return(results)
  }
}

#' Link RETSI and GETSI
#'
#' @param retsi RETSI calculated from DAR.
#' @param getsi GETSI calculated from DEG.
#' @inheritParams SPICEY
#' @return A \code{data.frame} combining RETSI and GETSI information
link_spicey <- function(retsi = NULL,
                        region_id = NULL,
                        getsi = NULL,
                        gene_id = NULL,
                        annotation = NULL) {
  
  #TODO: Check if this is working correctly!
  links <- retsi |> 
    dplyr::inner_join(annotation, by = c(region_id)) |> 
    dplyr::inner_join(getsi, by = c(gene_id, "cell_type"),
                      suffix = c("_ATAC", "_RNA"))
  
  return(links)
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
