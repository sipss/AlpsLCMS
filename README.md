# AlpsLCMS <img src='man/figures/AlpsLCMSlogo.png' align="right" height="139" />


[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Build Status](https://github.com/sipss/AlpsLCMS/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/sipss/AlpsLCMS/actions/) 
[![codecov.io](https://codecov.io/github/sipss/AlpsLCMS/coverage.svg?branch=master)](https://codecov.io/github/sipss/AlpsLCMS)
[![Documentation](https://img.shields.io/badge/documentation-pkgdown-informational)](https://sipss.github.io/AlpsLCMS/)

The goal of `AlpsLCMS` is to offer a data analysis preprocessing pipeline for LC/MS
metabolomic samples.

## Installation

AlpsLCMS can be installed with the `devtools` package. For this is needed Rtools and note that it uses packages from CRAN, from BioConductor and from git repositories: 

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
    
BiocManager::install(c("xcms", "mzR"))
BiocManager::install(c("CAMERA", "MAIT"))
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install(c("IPO", "faahKO"))
remotes::install_github("cbroeckl/RAMClustR")
devtools::install_git("https://gitlab.com/CarlBrunius/StatTools.git")

devtools::install_github("sipss/AlpsLCMS")
```


Quick start
========================

Checkout the [Introduction to AlpsLCMS](https://sipss.github.io/AlpsLCMS/articles/introduction-to-alpslcms.html) vignette that shows how to import data and preprocess it using `AlpsLCMS`.
