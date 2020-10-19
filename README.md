# NIHSlcms

<!-- badges: start -->
    [![Travis build status](https://travis-ci.com/Francisco-madrid-gambin/NIHSlcms.svg?branch=master)](https://travis-ci.com/Francisco-madrid-gambin/NIHSlcms)
<!-- badges: end -->
    <!-- badges: start -->
    [![Travis build status](https://travis-ci.com/Francisco-madrid-gambin/NIHSlcms.svg?branch=master)](https://travis-ci.com/Francisco-madrid-gambin/NIHSlcms)
<!-- badges: end -->
[![Build Status](https://github.com/sipss/AlpsNMR/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/sipss/AlpsNMR/actions/) [![codecov.io](https://codecov.io/github/sipss/AlpsNMR/coverage.svg?branch=master)](https://codecov.io/github/sipss/AlpsNMR) [![Documentation](https://img.shields.io/badge/documentation-pkgdown-informational)](https://sipss.github.io/AlpsNMR/)

The goal of `NIHSlcms` is to offer a data analysis preprocessing pipeline for LC/MS
metabolomic samples.

## Installation

NIHSlcms can be installed with the `devtools` package. For this is needed Rtools and note that it uses packages from CRAN, from BioConductor and from git repositories: 

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

)
BiocManager::install("xcms")
BiocManager::install("mzR")
BiocManager::install("CAMERA")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install("MAIT")
BiocManager::install("IPO")
devtools::install_github("cbroeckl/RAMClustR"
remotes::install_local("NIHSlcms_0.0.0.9007.tar.gz")

devtools::install_github("Francisco-madrid-gambin/NIHSlcms")
```


