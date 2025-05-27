library(dplyr)
library(tidyr)
library(GenomicRanges)
library(InteractionSet)
library(regioneR)
library(stringr)

# 1. Parse peak string ("chr-start-end") to GRanges
parse_peak <- function(peak_str) {
  parts <- str_split(peak_str, "-", simplify = TRUE)
  GRanges(seqnames = parts[, 1],
          ranges = IRanges(as.numeric(parts[, 2]),
                           as.numeric(parts[, 3])))
}

# 2. Build GInteractions object from links_df
make_links <- function(links_df, coaccess_threshold = 0.3) {
  links_df <- links_df %>%
    filter(coaccess > coaccess_threshold,
           !is.na(Peak1), !is.na(Peak2))

  anchor1 <- parse_peak(links_df$Peak1)
  anchor2 <- parse_peak(links_df$Peak2)

  links <- GInteractions(anchor1, anchor2)
  mcols(links)$coaccess <- links_df$coaccess
  mcols(links)$CCAN1 <- links_df$CCAN1
  mcols(links)$CCAN2 <- links_df$CCAN2
  links
}

# 3. Get CCAN membership for elements (REs or genes)
get_ccan <- function(links, gr, name_column, split = c("ccan", "name")) {
  split <- match.arg(split)
  hits1 <- findOverlaps(anchors(links, "first"), gr)
  hits2 <- findOverlaps(anchors(links, "second"), gr)

  ccan_df <- rbind(
    data.frame(ccan = mcols(links)$CCAN1[queryHits(hits1)],
               name = mcols(gr)[[name_column]][subjectHits(hits1)]),
    data.frame(ccan = mcols(links)$CCAN2[queryHits(hits2)],
               name = mcols(gr)[[name_column]][subjectHits(hits2)])
  ) |>
    unique() |>
    filter(!is.na(ccan), !is.na(name))

  if (split == "ccan") {
    split(ccan_df$name, ccan_df$ccan)
  } else {
    split(ccan_df$ccan, ccan_df$name)
  }
}

# 4. Main function to annotate RETSI peaks with coaccessible genes
annotate_with_coaccessibility <- function(links, retsi, getsi,
                                          name_column_peaks = "region",
                                          name_column_genes = "symbol") {
  # Step 1: Link REs to CCANs
  re_names <- mcols(retsi)[[name_column_peaks]]
  re_ccan <- get_ccan(links, retsi,
                      name_column = name_column_peaks, split = "name")
  re_ccan_vec <- re_ccan[re_names]
  re_ccan_vec[is.na(re_ccan_vec)] <- NA
  retsi$CCAN <- re_ccan_vec

  # Step 2: Link genes to CCANs
  gene_ccan <- get_ccan(links, getsi, name_column = name_column_genes, split = "ccan")
  used_ccans <- unique(unlist(re_ccan_vec, use.names = FALSE))
  gene_ccan <- gene_ccan[names(gene_ccan) %in% used_ccans]

  # Step 3: Add coaccessible genes to retsi
  retsi$genes_coacc <- lapply(retsi$CCAN, function(ccans) {
    if (length(ccans) == 0 || all(is.na(ccans))) return(NA_character_)
    unique(unlist(gene_ccan[as.character(ccans)], use.names = FALSE))
  })

  # Step 4: Expand each gene-coaccessible pair into separate GRanges
  keep <- lengths(retsi$genes_coacc) > 0
  gr_list <- rep(retsi[keep], lengths(retsi$genes_coacc[keep]))
  mcols(gr_list)$genes_coacc <- unlist(retsi$genes_coacc[keep], use.names = FALSE)
  return(gr_list)
}
