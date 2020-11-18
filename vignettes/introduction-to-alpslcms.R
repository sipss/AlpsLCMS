## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 7,
  fig.height = 5,
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(AlpsLCMS)

## ----message=FALSE, warning=FALSE---------------------------------------------
# Set the path of the example dataset. 
# Only 4 samples will be considered.
library(faahKO)
path <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
            recursive = TRUE)[5:8]

# Do not run if .mzXML format is used
polarity <- 1 # 1 for positive mode, 0 for negative mode

## ----Data wrangling-----------------------------------------------------------
# Set the mode to "onDisk" when you apply this function. Otherwise, it could
# take too much memory.
dataset <- lcms_read_samples(path, mode = "onDisk")

# Set the polarity
dataset@featureData@data[["polarity"]] <- rep(polarity, 
                                              length(dataset@featureData@data[["polarity"]]))
head(dataset)

## -----------------------------------------------------------------------------
# Create a metadata for this example dataset
metadata <- data.frame(sampleNames = basename(path),
                 treatment = c(rep("ko",2),
                               rep("wt",2)),
                 stringsAsFactors = FALSE)

## -----------------------------------------------------------------------------
dataset <- lcms_meta_add(dataset, metadata, by = "sampleNames")
phData(dataset)

## -----------------------------------------------------------------------------
tics <- lcms_tics(dataset, treatment = "treatment")

lcms_plot_tics(tics,
               treatment = treatment, 
               plot_type = "spec")

lcms_plot_tics(tics, treatment = treatment,
               plot_type = "boxplot")

## -----------------------------------------------------------------------------
# Range of the retention time (minutes) to include in further analyses
rt = c(50, 60)
ms = c(200, 500)

## -----------------------------------------------------------------------------
dataset_shorter <- lcms_filter_rt_min(dataset, rt = rt)
dataset_shorter <- lcms_filter_mz(dataset_shorter, mz = ms)
tics <- lcms_tics(dataset_shorter, treatment = "treatment")

lcms_plot_tics(tics,
               treatment = treatment, 
               plot_type = "spec")

## ----Creating parameters------------------------------------------------------
# if to optimize or not
optimize <- TRUE
nSlaves <- 1
subset <- dataset[3:4]
classes <- dataset@phenoData@data[["treatment"]]

