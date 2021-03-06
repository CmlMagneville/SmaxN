
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SmaxN

The `SmaxN` package provides a unique function to compute maximisation
of abundances index. It notably helps to compute the Synchronised maxN
(SmaxN) based on the use of temporally synchronised cameras.

## Installation

You can install the `SmaxN` package from [GitHub](https://github.com/)
with:

``` r
## Install < remotes > package (if not already installed) ----
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

## Install < SmaxN > package from GitHub ----
remotes::install_github("CmlMagneville/SmaxN", build_vignettes = TRUE)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SmaxN)

# Build distance dataframe for the example:
dist_df_ex <- data.frame("A" = c(0, 2, 5, 5), "B" = c(2, 0, 5, 5), 
                         "C" = c(5, 5, 0, 4), "D" = c(5, 5, 4, 0))
rownames(dist_df_ex) <- c("A", "B", "C", "D")

# Build distance dataframe for the example:
abund_df_ex <- data.frame("A" = c(0, 1, 3, 7, 2, 2, 3, 0, 6, 2, 0, 1), 
                          "B" = c(2, 2, 2, 2, 0, 0, 0, 0, 1, 2, 1, 0), 
                          "C" = c(2, 0, 1, 0, 0, 4, 2, 2, 3, 0, 0, 4), 
                          "D" = c(0, 1, 0, 1, 0, 6, 1, 1, 6, 4, 2, 1))
 
# Run the general function of the package :
SmaxN::compute.max.abund(dist_df = dist_df_ex, 
                 fish_speed = 1.6, 
                 abund_df   = abund_df_ex)
#> $SmaxN
#> [1] 19
#> 
#> $maxN
#> [1] 7
#> 
#> $maxN_cam
#> A B C D 
#> 7 2 4 6 
#> 
#> $maxN_row
#>  [1] 2 2 3 7 2 6 3 2 6 4 2 4
#> 
#> $SmaxN_row
#> [1] 16
```
