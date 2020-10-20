# NIHSlcms

[![Build Status](https://github.com/sipss/NIHSlcms/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/sipss/NIHSlcms/actions/) 
[![codecov.io](https://codecov.io/github/sipss/NIHSlcms/coverage.svg?branch=master)](https://codecov.io/github/sipss/NIHSlcms)
[![Documentation](https://img.shields.io/badge/documentation-pkgdown-informational)](https://sipss.github.io/NIHSlcms/)

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
Quick start
=============

Checkout the [Introduction to NIHSlcms](https://sipss.github.io/NIHSlcms/articles/introduction-to-alpsnmr.html) vignette that shows how to import data and preprocess it using NIHSlcms.

See also the [tutorial](https://github.com/sipss/AlpsNMR/blob/master/vignettes/tutorial.pdf) with a real dataset from beginning to end, including all the steps of untargeted metabolomics analysis. To run the [tutorial](https://github.com/sipss/AlpsNMR/blob/master/vignettes/tutorial.pdf), you can download the MTBLS242 dataset from the public [MetaboLights repository](https://www.ebi.ac.uk/metabolights/MTBLS242), or download and unzip the contents (spectra and metadata) of this [Dropbox link](https://dl.dropboxusercontent.com/s/0snivrsd7m82yey/MTBLS242.zip?dl=0).

