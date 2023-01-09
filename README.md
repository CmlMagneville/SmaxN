
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

In this example, we see how to use the `SmaxN.computation` function of
the `SmaxN` package which computes the SmaxN metric and the maxN metrics
among other index (for a better explanation, don’t hesitate to look at
the [“How to compute the SmaxN and other abundance metrics?
tutorial”]()) on the [SmaxN
website](https://cmlmagneville.github.io/SmaxN/).

``` r
# Build distance dataframe for the example:
dist_df_ex <- data.frame("A" = c(0, 2, 5, 5), "B" = c(2, 0, 5, 5), 
                         "C" = c(5, 5, 0, 4), "D" = c(5, 5, 4, 0))
rownames(dist_df_ex) <- c("A", "B", "C", "D")

# Build distance dataframe for the example:
abund_df_ex <- data.frame("A" = c(0, 1, 3, 7, 2, 2, 3, 0, 6, 2, 0, 1), 
                          "B" = c(2, 2, 2, 2, 0, 0, 0, 0, 1, 2, 1, 0), 
                          "C" = c(2, 0, 1, 0, 0, 4, 2, 2, 3, 0, 0, 4), 
                          "D" = c(0, 1, 0, 1, 0, 6, 1, 1, 6, 4, 2, 1))
# Call the package:
library("SmaxN")
 
# Run the general function of the package :
SmaxN_results <- SmaxN::SmaxN.computation(abund_df = abund_df_ex, 
                                          speed = 1.6, 
                                          dist_df = dist_df_ex)
#> [1] "!!!!!! Starting row  1 on 2"
#> [1] "50%"
SmaxN_results
#> $maxN
#> [1] 7
#> 
#> $SmaxN
#> [1] 19
#> 
#> $path_saved
#> $path_saved[[1]]
#> $path_saved[[1]][[1]]
#>   value cam_nm timestep
#> 1     7      A        4
#> 2     2      B        4
#> 3     4      C        6
#> 4     6      D        6
#> 
#> 
#> 
#> $number_SmaxN_path
#> [1] 1
#> 
#> $SmaxN_timestep
#> [1] 16
#> 
#> $maxN_cam
#> A B C D 
#> 7 2 4 6 
#> 
#> $maxN_timestep
#>  1  2  3  4  5  6  7  8  9 10 11 12 
#>  2  2  3  7  2  6  3  2  6  4  2  4
```

## Citation

Please cite this package as:

*Magneville et al.* (2022). SmaxN: Maximisation of abundances using
synchronised cameras in R. R package version 0.2.
<https://github.com/CmlMagneville/SmaxN>
