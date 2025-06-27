# SPICEY <img src="man/figures/logo_spicey.png" width="140px" height="140px" align="right" style="padding-left:10px;background-color:white;"/>

<!-- badges: start -->

<!-- badges: end -->

## Bioconductor release status

| Branch | R CMD check | Last updated |
|:------------------------:|:------------------------:|:-------------------:|
| [*devel*](http://bioconductor.org/packages/devel/bioc/html/SPICEY.html) | [![Bioconductor-devel Build Status](http://bioconductor.org/shields/build/devel/bioc/SPICEY.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/SPICEY) | ![](http://bioconductor.org/shields/lastcommit/devel/bioc/SPICEY.svg) |
| [*release*](http://bioconductor.org/packages/release/bioc/html/SPICEY.html) | [![Bioconductor-release Build Status](http://bioconductor.org/shields/build/release/bioc/SPICEY.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/SPICEY) | ![](http://bioconductor.org/shields/lastcommit/release/bioc/SPICEY.svg) |

The goal of SPICEY is to provide a user-friendly pipeline for quantifying and visualizing tissue specificity specificity

## Installation

You can install the latest release of `SPICEY` the github repository:

```         
devtools::install_github("georginafp/SPICEY")
```

Now you can load the package using `library(SPICEY)`.

## Basic usage

For detailed instructions on how to use SPICEY, please see the [vignette](https://georginafp.github.io/SPICEY/articles/SPICEY.html).

``` r
library(SPICEY)
```

``` r
# Compute GETSI 
results <- run_spicey(rna = rna)

# Compute RETSI
results <- run_spicey(atac=atac)

# Compute GETSI + RETSI
results <- run_spicey(atac=atac, rna=rna)

# Compute GETSI + RETSI and link RE to target genes through nearest gene method
results <- run_spicey(rna = rna, atac=atac, link_method = "nearest")

# Compute GETSI + RETSI and link RE to target genes through coaccessibility method
results <- run_spicey(
  atac = atac_data,
  rna = rna_data,
  link_method = "coaccessibility",
  links = coaccessibility_links,
  coaccess_cutoff_override = 0.25,  # default
  filter_promoter_distal = TRUE     # default
)
```

## Code of Conduct

Please note that the SPICEY project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
