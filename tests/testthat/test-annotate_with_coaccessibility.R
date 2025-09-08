data(atac)
data(cicero_links)

test_that("annotate_with_coaccessibility returns GRanges with gene links", {
  peaks <- unique(unlist(atac)[, "region_id"])
  ann <- annotate_with_coaccessibility(
    peaks = peaks,
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    links_df = cicero_links,
    annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
    protein_coding_only = TRUE,
    verbose = TRUE,
    add_tss_annotation = FALSE,
    upstream = 2000,
    downstream = 2000)
  expect_s3_class(ann, "data.frame")
  expect_true(all(c("region_id", "gene_id") %in% colnames(ann)))
})
