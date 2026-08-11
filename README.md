# dynclusts <img src="man/figures/dynclusts.png" alt="dynclusts" width="140" height="150" align="right"> 

<!-- badges: start -->

[![R-CMD-check](https://github.com/marinsantiago/dynclusts/workflows/R-CMD-check/badge.svg)](https://github.com/marinsantiago/dynclusts/workflows/R-CMD-check/badge.svg)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

</br>

## Overview

The R package `dynclusts` (developer's version) performs dynamic clustering through an 
[autoregressive logistic-beta Stirling-gamma process](https://arxiv.org/abs/2601.04625) 
as described in Marin et al. (2026+).

## Installation

You can install the latest developer's version via `pak` as:

``` r
# install.packages("pak")
pak::pak("marinsantiago/dynclusts")
```

On the other hand, if you wish to install the package from the `dynclusts.zip` file in the supplementary materials to Marin et al. (2026+):

  1. In R, set your working directory to the folder `dynclusts`.
  
  2. Run the following R code:
  
``` r
# install.packages("devtools")
devtools::build()
devtools::install()
```

## Usage

Detailed guidelines for using the package functions are referred to their help pages in R. Additional examples are available at [https://github.com/marinsantiago/dynclusts-applications](https://github.com/marinsantiago/dynclusts-applications).

## <a name="cite"></a> Citation

If you use any part of this code in your work, please consider citing our paper:

```
@misc{marin_dynclusts,
  title         = {{B}ayesian nonparametric modeling of dynamic pollution clusters through an autoregressive logistic-beta {S}tirling-gamma process}, 
  author        = {Santiago Marin and Bronwyn Loong and Anton H. Westveld},
  year          = {2026},
  eprint        = {2601.04625},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ME}
}
```

## <a name="refs"></a> References

Marin, S., Loong, B., and Westveld, A. H. (2026+), "Bayesian nonparametric modeling of dynamic pollution clusters through an autoregressive logistic-beta Stirling-gamma process"
