CHANGES IN VERSION 0.99.5
--------------------------

  o Updated README file with new example of heatmap plot.

CHANGES IN VERSION 0.99.4
--------------------------

  o Updated spicey_plot.R: the spicey_heatmap() function no longer computes z-scores from SPICEY values; it now directly visualizes the original SPICEY scores.
  o Updated the vignette accordingly to remove all references and calculations related to z-scores.

CHANGES IN VERSION 0.99.3
--------------------------

  o Minor update: modified one line in the internal plot_heatmap() function to adjust the color gradient for visualization purposes. 

CHANGES IN VERSION 0.99.2
--------------------------

  o Minor update: Fix the doi issue for citation.
  
CHANGES IN VERSION 0.99.1
--------------------------

  o Updated R dependency (R ≥ 4.5)
  o Added unit tests for main functions following Bioconductor guidelines.
  o Added fnd role in Authors@R (FPI fellowship acknowledgment).
  o Verified all :: usages; kept for clarity with proper Imports: tracking.
  o Refactored redundant code into helper functions (.standardize_peaks, .add_tss_annotation); streamlined heatmap logic.
  o Added package man page.
  o Updated vignette: removed GitHub installation, added Bioconductor installation via BiocManager.
  o Fixed minor typo in specificity_index.Rd.
  
CHANGES IN VERSION 0.99.0
--------------------------

  o First public release of SPICEY.
  o Added a `NEWS.md` file to track changes to the package.
