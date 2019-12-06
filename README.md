
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GSVAutils

<!-- badges: start -->

<!-- badges: end -->

The goal of GSVAutils is to provide severall functions used in TIGS
project.

## Installation

``` r
remotes::install_github("GSVAutils")
```

## Example

This is a basic example which shows you how to run GSVA.

``` r
library(GSVAutils)
#> Warning: replacing previous import 'Biobase::combine' by 'dplyr::combine' when
#> loading 'GSVAutils'
## basic example code
data("example_gsets")
data("example_data")

ExprList = GenTibbleList(list(example_data))
#> Find nothing about gene symbol in fData, try search it...
#> Find gene_assignment
#> Processing...
#> Done.
res <- ApplyGSVA(example_gsets,
  group_col = "Cell_type",
  gene_col = "Symbol", ExprMatList = ExprList, method = "gsva")
#> Estimating GSVA scores for 21 gene sets.
#> Computing observed enrichment scores
#> Estimating ECDFs with Gaussian kernels
#> Using parallel with 8 cores
#> 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |======================================================================| 100%
```

More about GSVA please run `?GSVA::gsva` in R console.
