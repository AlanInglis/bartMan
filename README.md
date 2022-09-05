
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bartMan

bartMan is an R-package for investigating and visualising Bayesian
Additive Regression Tree (BART) model fits. We construct conventional
plots to analyze a model’s performance and stability as well as create
new tree-based plots to analyze variable importance, interaction, and
tree structure. We employ Value Suppressing Uncertainty Palettes (VSUP)
to construct heatmaps that display variable importance and interactions
jointly using color scale to represent posterior uncertainty. Our
visualizations are designed to work with the most popular BART R
packages available, namely BART, dbarts, and bartMachine. A practical
example of the package in use can be found here:
<https://alaninglis.github.io/bartMan/articles/bartManVignette.html>

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
