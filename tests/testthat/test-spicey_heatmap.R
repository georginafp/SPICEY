data(rna)
data(atac)
data(cicero_links)

test_that("spicey_heatmap runs for RETSI, GETSI, and SPICEY", {
  peaks <- unique(unlist(atac)[, c("region_id")])
  annotation_coacc <- annotate_with_coaccessibility(
    peaks = peaks,
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    links_df = cicero_links,
    annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
    protein_coding_only = TRUE,
    verbose = TRUE,
    add_tss_annotation = FALSE,
    upstream = 2000,
    downstream = 2000)

  spicey_coacc <- SPICEY(
    rna = rna,
    atac = atac,
    annotation = annotation_coacc)

  spicey_coacc$RETSI <- spicey_coacc$RETSI |> dplyr::left_join(annotation_coacc, by = c("region_id"))

  p1 <- spicey_heatmap(spicey_coacc$RETSI, spicey_measure = "RETSI")
  expect_s3_class(p1, "ggplot")

  p2 <- spicey_heatmap(spicey_coacc$GETSI, spicey_measure = "GETSI")
  expect_s3_class(p2, "ggplot")

  p3 <- spicey_heatmap(spicey_coacc$linked, spicey_measure = "SPICEY", combined_score = FALSE)
  expect_true(inherits(p3, "patchwork") || inherits(p3, "gg"))

  p4 <- spicey_heatmap(spicey_coacc$linked, spicey_measure = "SPICEY", combined_score = TRUE)
  expect_s3_class(p4, "ggplot")
})
