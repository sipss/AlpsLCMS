# NIHSlcms

The goal of NIHSlcms is to offer a data analysis preprocessing pipeline for LC/MS
metabolomic samples.

## Installation

You can install the released version of NIHSlcms with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

devtools::install_github("sneumann/mzR")
BiocManager::install("xcms")
BiocManager::install("CAMERA")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install("MAIT")
BiocManager::install("IPO")
remotes::install_local("NIHSlcms_0.0.0.9007.tar.gz")
```


