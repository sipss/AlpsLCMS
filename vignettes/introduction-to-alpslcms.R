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
library(faahKO)
path <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
            recursive = TRUE)
polarity <- 1 # 1 for positive mode, 0 for negative mode

## ----Data wrangling-----------------------------------------------------------
# Be careful setting the mode to "onDisk" when you apply this function.
dataset <- lcms_read_samples(path, mode = "onDisk")
dataset@featureData@data[["polarity"]] <- rep(polarity, length(dataset@featureData@data[["polarity"]]))
head(dataset)

## -----------------------------------------------------------------------------
metadata <- data.frame(sampleNames = basename(path),
                 treatment = c(rep("ko",6),
                               rep("wt",6)),
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
prep_parm_p <- NULL
prep_parm_p$peakwidth <- c(20, 80)
prep_parm_p$noise <- 5000
prep_parm_p$snthresh <- 3
prep_parm_p$prefilter <- c(6, 5000)
prep_parm_p$centerSample <- "wMean"
prep_parm_p$integrate <- 2
prep_parm_p$mzdiff <- -0.001
prep_parm_p$profStep <- 0.005
prep_parm_p$minFraction <- 0.2
prep_parm_p$ppm <- 25
prep_parm_p$mzCenterFun <- "wMean"
prep_parm_p$fitgauss <- FALSE
prep_parm_p$verbose.columns <- FALSE

classes <- dataset@phenoData@data[["treatment"]]

## ----Peak detection-----------------------------------------------------------
peakdet = find_peaks_cwp(dataset, 
                         params = prep_parm_p)

message("Number of detected peaks")
peakdet@msFeatureData[["chromPeakData"]]@nrows
message("")

message("Parameters")
peakdet@.processHistory[[1]]@param

xcms::plotChromPeakImage(peakdet)

## ----Correspondence-----------------------------------------------------------
new_params <- PeakDensityPar(sampleGroups = classes, 
                               binSize = 0.6)

peakgrouped <- groupPeaks(peakdet, 
                         param = new_params)

## ----Alignment----------------------------------------------------------------
pgp <- PeakGroupsPar(minFraction = 0.8,
                       extraPeaks = 1, 
                       smooth = "loess",
                       span = 0.4,
                       family = "gaussian")

## Get the peak groups that would be used for alignment.
xdata_aling <- adjustRT(peakgrouped, param = pgp)

rt_plot = lcms_retention_time_alignment_plot(xdata_aling)
rt_plot

## REGROUPING
new_params <- PeakDensityPar(sampleGroups = classes,
                             bw = 30, 
                             minFraction = 0.4)

peakgrouped = groupPeaks(xdata_aling, param = new_params)

## -----------------------------------------------------------------------------
lcms_plot_chrom_peak_image(peakdet, binSize = 5,
                           xlim = NULL,
                           log = FALSE,
                           xlab = "retention time (min)",
                           yaxt = par("yaxt"),
                           main = "Detected Peaks (unprocessed)")

lcms_plot_chrom_peak_image(peakgrouped, binSize = 5,
                           xlim = NULL,
                           log = FALSE,
                           xlab = "retention time (min)",
                           yaxt = par("yaxt"),
                           main = "Detected Peaks (processed)")

## ----Imputation I-------------------------------------------------------------
message("Missing values found in the processed dataset: ", sum(is.na(featureValues(peakgrouped))))

peakgrouped_imp <- lcms_fill_chrom_peaks(peakgrouped)
cat("Imputing values...\n")

message("Missing values found after fill_chrom_peaks: ", sum(is.na(featureValues(peakgrouped_imp))))

## ----Feature table------------------------------------------------------------
xdata = featureValues(peakgrouped_imp,
                      method = "maxint",
                      value = "into",
                      filled = TRUE, 
                      missing = "rowmin_half")
xdata <-  t(xdata)
feature <- featureDefinitions(peakgrouped_imp)
feature <- feature@listData
featNames <- paste0(feature$mzmed,"_",feature$rtmed)
colnames(xdata) <- featNames

message("Missing values in the feature table: ",
sum(is.na(xdata)))

## ----echo = FALSE-------------------------------------------------------------
# xdataImp <- xdata
# xdataImputed <- as.data.frame(xdataImp, stringsAsFactors = FALSE)

# Get mz and rt columns for the feature table
mz <- colnames(xdata) %>%
  stringr::str_split(.,"\\_") %>% 
  lapply(.,function(x) x[1]) %>% 
  unlist() %>% 
  as.numeric()

#You can get rt also
rt <- colnames(xdata) %>%
  stringr::str_split(.,"\\_") %>% 
  lapply(.,function(x) x[2]) %>% 
  unlist() %>% 
  as.numeric()
rt <- rt/60

## ----Params Data reduction----------------------------------------------------
st <- getRamSt(peakgrouped_imp)
sr <- 0.6

#List of adducts for do.findmain
#adducts_list = c("[M+H-H2O]+")
adducts_list = c()

## Building the defineExperiment manually
## Change for your convenience (e.g. GC-MS)
value <- c(rep("fill", 4), "LC-MS")
design <- as.data.frame(value)
rownme <- c("Experiment", "Species", "Sample",
            "Contributer", "platform")  
rownames(design) <- rownme

value <- c(rep("fill", 13), "1")
instrument <- as.data.frame(value)
rownm <- c("chrominst", "msinst", "column", 
           "solvA", "solvB", "CE1", "CE2", 
           "mstype", "msmode", "ionization", 
           "colgas", "msscanrange", "conevol", 
           "MSlevs")
rownames(instrument) <- rownm

Experiment <- list(design =  design, instrument = instrument)

## ----warning=FALSE------------------------------------------------------------
RC <- clustering(xcmsObj = peakgrouped_imp,
                 featdelim = ".",
                 st = st,
                 sr = sr,
                 ExpDes = Experiment,
                 normalize = "TIC",
                 deepSplit = TRUE,
                 sampNameCol = 1,
                 mspout = FALSE,
                 fftempdir = getwd())

RC <- do.findmain(RC,
                  nls = adducts_list,
                  mode = "positive",
                  mzabs.error = 0.005,
                  ppm.error = 5,
                  plot.findmain = FALSE,
                  writeMat = FALSE,
                  writeMS = FALSE)


## ----Reduced feature table----------------------------------------------------
labeled_adducts <- labelled_adducts(RC)
representative_ions <- labeled_adducts$Representative_ions
xdata_reduced <- feature_reduction(xdata, representative_ions, RC)

## -----------------------------------------------------------------------------
stat <- function(x){stats::wilcox.test(x ~ classes, xdata_reduced)$p.value}
abcd <- data.frame(apply(FUN = stat,
                         MARGIN = 2,
                         X = xdata_reduced))
colnames(abcd) <- c("p_Wilc")
abcd$id <- row.names(abcd)
result <- abcd
summary(result$p_Wilc)
message("\nNumber of features < 0.05 nominal p-value ", 
sum(result$p_Wilc < 0.05))
head(result[result$p_Wilc < 0.05,1])

#FDR
fdr.wilcox <- stats::p.adjust(result$p_Wilc, method = "fdr")
result <- cbind(result, fdr.wilcox)
message("\nNumber of features fdr-corrected p value of < 0.05 is ", 
sum(result$fdr.wilcox < 0.05))

## ----Annotation univariate postivie-------------------------------------------
# We use the selected vips
univ_feat <- result[result$fdr.wilcox < 0.05,"id"]

# Untargeted assignation
# Creating mz column
mzr <- univ_feat %>%
  stringr::str_split(.,"\\_") %>%
  lapply(.,function(x) x[1]) %>%
  unlist() %>%
  as.numeric()

all.equal(length(univ_feat), length(mzr))

tdata_reduced_univ <- data.frame(mz = mzr)
if(dim(tdata_reduced_univ)[1]>0){
  result_POS_HMDB_univ <- assignation_pos_HMDB(tdata_reduced_univ)
  head(result_POS_HMDB_univ)
} else {
  message("There is no significant features in Wilcox test")
}

## ----Annotation of all features positive--------------------------------------
# Untargeted assignation
# Creating mz column

mzr <- colnames(as.matrix(xdata)) %>%
  stringr::str_split(.,"\\_") %>%
  lapply(.,function(x) x[1]) %>%
  unlist() %>%
  as.numeric()

all.equal(dim(xdata)[2], length(mzr))

tdata_reduced_all_features <- data.frame(mz = mzr)
result_POS_HMDB_all_features <- assignation_pos_HMDB(tdata_reduced_all_features)
head(result_POS_HMDB_all_features)

