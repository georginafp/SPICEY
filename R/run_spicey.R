run_spicey <- function(atac = NULL,
                       rna = NULL,
                       annot_method = NULL,
                       links = NULL,
                       link_spicey_measures = FALSE,
                       coaccess_cutoff_override = 0.25,
                       filter_promoter_distal = TRUE,
                       filter_protein_coding = TRUE,
                       txdb = NULL,
                       keep_mito = FALSE,
                       annot_dbi = NULL,
                       add_tss_annotation = FALSE,
                       verbose = TRUE) {

  if (is.null(atac) && is.null(rna)) {
    stop("At least one of 'atac' or 'rna' data must be provided.")
  }

  # Compute RETSI
  if (!is.null(atac)) {
    retsi <- spicey_retsi(atac)
  }

  # Compute GETSI
  if (!is.null(rna)) {
    getsi <- spicey_getsi(rna)
  }

  # If no annotation method, return based on inputs
  if (is.null(annot_method)) {
    if (!is.null(atac) && is.null(rna)) {
      message("🌶 SPICEY pipeline successfully completed")
      return(retsi)
    }
    if (is.null(atac) && !is.null(rna)) {
      message("🌶 SPICEY pipeline successfully completed")
      return(getsi)
    }
    if (!is.null(atac) && !is.null(rna)) {
      message("🌶 SPICEY pipeline successfully completed")
      return(list(retsi = retsi, getsi = getsi))
    }
  }

  # Annotation requires ATAC data
  if (is.null(atac)) {
    stop("ATAC data must be provided for region-to-gene annotation.")
  }

  if (!(annot_method %in% c("nearest", "coaccessibility"))) {
    stop("annot_method must be either 'nearest' or 'coaccessibility'.")
  }

  if (annot_method == "nearest") {
    retsi_annotated <- annotate_with_nearest(
      retsi = retsi,
      txdb = txdb,
      keep_mito = keep_mito,
      annot_dbi = annot_dbi,
      protein_coding_only = filter_protein_coding,
      verbose = verbose,
      add_tss_annotation = add_tss_annotation
    )
  }

  if (annot_method == "coaccessibility") {
    if (is.null(links)) {
      stop("Links must be provided when using coaccessibility annotation method.")
    }
    retsi_annotated <- annotate_with_coaccessibility(
      links = links,
      retsi = retsi,
      txdb = txdb,
      keep_mito = keep_mito,
      annot_dbi = annot_dbi,
      protein_coding_only = filter_protein_coding,
      verbose = verbose,
      coaccess_cutoff_override = coaccess_cutoff_override,
      filter_promoter_distal = filter_promoter_distal,
      add_tss_annotation = add_tss_annotation
    )
  }

  # Now handle linking or returning annotated + getsi

  if (link_spicey_measures) {
    if (is.null(rna)) {
      stop("RNA data must be provided for linking RETSI to gene expression (GETSI).")
    }
    if (verbose) message("→ Computing GETSI (if not already computed)...")
    if (!exists("getsi")) {
      getsi <- spicey_getsi(rna)
    }

    if (verbose) message("→ Linking RETSI and GETSI via ", annot_method, "...")

    combined <- NULL
    if (annot_method == "nearest") {
      combined <- link_spicey_nearest(
        retsi_annotated_nearest = retsi_annotated,
        getsi = getsi
      )
    } else {
      combined <- link_spicey_coaccessible(
        retsi_annotated_coacc = retsi_annotated,
        getsi = getsi
      )
    }
    # Return the combined linked data.frame/tibble, NOT a list
    message("🌶 SPICEY pipeline successfully completed")
    return(combined)
  } else {
    # Linking is FALSE: return a list with retsi_annotated and getsi if both available
    if (!is.null(rna)) {
      message("🌶 SPICEY pipeline successfully completed")
      return(list(retsi_annotated = retsi_annotated, getsi = getsi))
    } else {
      message("🌶 SPICEY pipeline successfully completed")
      return(retsi_annotated)
    }
  }
}
