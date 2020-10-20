## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 7,
  fig.height = 5,
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(NIHSlcms)
library(ggplot2)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(faahKO)
path <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
            recursive = TRUE)[c(1:3, 7:9)]
polarity <- 1 # 1 for positive mode, 0 for negative mode

## ----Data wrangling-----------------------------------------------------------
# Be careful setting the mode to "onDisk" when you apply this function.
dataset <- readMSData(path, mode = "onDisk")
dataset@featureData@data[["polarity"]] <- rep(polarity, length(dataset@featureData@data[["polarity"]]))
message("Is the polarity correctly loaded? \n")
cat("\n")
print(dataset)

## -----------------------------------------------------------------------------
metadata <- data.frame(sampleNames = basename(path),
                 treatment = c("ko", "ko", "ko", "wt", "wt", "wt"),
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

## -----------------------------------------------------------------------------
metabolite <- "C22:0"
retention_time <- NULL
mz <- 446.3309
dataset <- dataset

## -----------------------------------------------------------------------------
# retention_time_max <- (retention_time*60)+120
# retention_time_min <- (retention_time*60)-120
# rtr <- c(retention_time_min, retention_time_max)

mz_max <- mz + 0.01
mz_min <- mz - 0.01
mzr <- c(mz_min, mz_max)

## extract the chromatogram
chr_raw <- chromatogram(dataset,
                        # rt = rtr,
                        mz = mzr)

plot(chr_raw, col = "blue", main = metabolite)
# plot(chr_raw_leucines, col = "blue", xlim = c(330,370), ylim = c(0,1e+8), main = metabolite)
# plot(chr_raw_leucines[[2]], col = "blue", xlim = c(330,370), ylim = c(0,1e+8), main = metabolite)

## ----Optimize-----------------------------------------------------------------
opt_path <- tempdir(system.file("opt_path", package = "NIHSlcms"))
optimize <- TRUE
nSlaves <- 1

## -----------------------------------------------------------------------------
default_peakpeaking_params <- lcms_default_peakpicking_params(noise = c(5e+05, 1e+06),
                                                              snthresh = 3, 
                                                              min_peakwidth = c(5, 20), 
                                                              max_peakwidth = c(35, 60), 
                                                              optimize = optimize)
# 
# default_peakpeaking_params$prefilter <- 3 
# default_peakpeaking_params$value_of_prefilter <- 1000
# default_peakpeaking_params$integrate <- 1
# 
# resultPeakpicking <- lcms_peakpicking_optimization(
#   dataset = dataset,
#   peakpickingParameters = default_peakpeaking_params,
#   nSlaves = nSlaves,
#   opt_path = opt_path,
#   subdir = NULL,
#   plots = FALSE)

## -----------------------------------------------------------------------------
# optimizedXcmsSet <- resultPeakpicking$best_settings$xset


