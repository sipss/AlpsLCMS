# NIHSlcms

The goal of NIHSlcms is to offer a data analysis preprocessing pipeline for LC/MS
metabolomic samples.

## Installation

You can install the released version of NIHSlcms with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("Biobase")    
BiocManager::install("CAMERA")
BiocManager::install("MSnbase")
BiocManager::install("xcms")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install("MAIT")
BiocManager::install("IPO")
remotes::install_local("NIHSlcms_0.0.0.9004.tar.gz")
```

