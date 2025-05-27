library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(Seurat)
library(tidyr)
library(regioneR)
library(purrr)
library(Signac)
library(purrr)
library(ggside)
library(viridis)
library(tibble)
library(ggsignif)
library(RColorBrewer)
library(ggplotify)
library(pheatmap)
library(patchwork)
library(scales)
library(data.table)
library(here)
library(InteractionSet)
library(stringr)
library(tools)


source("R/getsi.R")
source("R/retsi.R")
source("R/utils.R")
source("R/link_spicey.R")


run_spicey <- function(atac_path, rna_path, links_path) {

  message("→ Step 1: Reading input data...")
  atac <- read_input_file(atac_path)
  rna <- read_input_file(rna_path)
  links_df <- read_input_file(links_path)

  message("→ Step 2: Creating GInteractions object from links...")
  links <- make_links(links_df, coaccess_threshold = 0.3)

  message("→ Step 3: Computing RETSI and RETSI entropy...")
  atac_scored <- compute_retsi(atac)
  atac_entropy <- compute_entropy_retsi(atac_scored)
  mcols(atac_scored)$region <- paste0(seqnames(atac_scored), ":", start(atac_scored), "-", end(atac_scored))
  atac_scored_df <- as.data.frame(atac_scored) %>%
    left_join(atac_entropy, by = "region") %>%
    select(-c("width", "strand")) %>%
    toGRanges()

  message("→ Step 4: Computing GETSI and GETSI entropy...")
  rna_scored <- compute_getsi(rna)
  rna_entropy <- compute_entropy_getsi(rna_scored)
  rna_scored_df <- as.data.frame(rna_scored) %>%
    left_join(rna_entropy, by = "symbol") %>%
    select(-c("width", "strand")) %>%
    toGRanges()

  message("→ Step 5: Linking REs and genes via co-accessibility...")
  gr_links <- annotate_with_coaccessibility(
    links = links,
    retsi = atac_scored_df,
    getsi = rna_scored_df,
    name_column_peaks = "region",
    name_column_genes = "symbol"
  )

  message("→ Step 6: Merging with GETSI scores...")
  getsi_df <- as.data.frame(rna_scored_df) %>%
    select(symbol, GETSI, cell_type, norm_entropy) %>%
    rename(GETSI_coacc = GETSI, GETSI_entropy = norm_entropy)

  final_df <- gr_links %>%
    data.frame() %>%
    rename(RETSI_entropy = norm_entropy) %>%
    left_join(getsi_df, by = c("genes_coacc" = "symbol", "cell_type")) %>%
    select(-c("avg_FC", "max_log2FC", "weight", "norm_log2FC", "entropy",
              "CCAN"))

  message("✅ SPICEY pipeline completed successfully.")
  return(final_df)
}


final_result <- run_spicey(
  atac_path = "data/FINAL_ATAC.rds",
  rna_path = "data/FINAL_RNA.rds",
  links_path = "data/COACC_LINKS.rds"
  )
saveRDS(final_result, "data/SPICEY_FINAL.rds")
