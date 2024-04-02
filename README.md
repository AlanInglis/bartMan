
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bartMan

<!-- badges: start -->

<img src="https://raw.githubusercontent.com/AlanInglis/bartMan/master/badge/bartmanLogo1.png" width="240" height="276" align="right" />
<!-- badges: end --> For more detailed information and a comprehensive
discussion, please refer to our paper associated with this document,
available at DOI: (<https://doi.org/10.52933/jdssv.v4i1.79>). bartMan is
an R-package for investigating and visualising Bayesian Additive
Regression Tree (BART) model fits. We construct conventional plots to
analyze a model’s performance and stability as well as create new
tree-based plots to analyze variable importance, interaction, and tree
structure. We employ Value Suppressing Uncertainty Palettes (VSUP) to
construct heatmaps that display variable importance and interactions
jointly using color scale to represent posterior uncertainty. Our
visualizations are designed to work with the most popular BART R
packages available, namely BART, dbarts, and bartMachine. A practical
example of the package in use can be found here:
(<https://alaninglis.github.io/bartMan/articles/bartManVignette.html>)

## Installation

You can install the development version from
[GitHub](https://github.com/AlanInglis/bartMan) with:

``` r
# install.packages("devtools")
devtools::install_github("AlanInglis/bartMan")
```

You can then load the package with:

``` r
library(bartMan)
```
