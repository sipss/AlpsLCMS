#' Get st parameter value for clustering
#'
#' getRamSt calculates the optimal value for the [AlpsLCMS::clustering]
#' function. It uses the `lcms_dataset` after xcms preprocessing.
#'
#' @param lcms_dataset A lcms_dataset
#' @return recommended st parameter value for [AlpsLCMS::clustering]
#' @export
getRamSt <- function(lcms_dataset) {
  featInfo <- xcms::featureDefinitions(lcms_dataset)
  histo <- graphics::hist((featInfo$rtmax-featInfo$rtmin)/2)
  st <- round(stats::median(featInfo$rtmax-featInfo$rtmin)/2,
              digits = 2)
  if(st == 0 ){
    st <- round(mean(featInfo$rtmax-featInfo$rtmin)/2,
                digits = 2)
  }
  if(st == 0){
    st = 0.01
  }
  graphics::abline(v=st)
  return(st)
}

#' Clustering
#'
#' `clustering` is a wrapper of the [RAMClustR::ramclustR] from `RAMClustR`
#' package. It performs a clustering of features with a given sigma for
#' retention time similarity `st` and for correlation similarity `sr`. Note
#' that, in addition to the `sr`, the argument `deepSplit = TRUE` might be
#' critical to avoid several metabolites in a single cluster.
#'
#' @inheritParams RAMClustR::ramclustR
#' @inherit RAMClustR::ramclustR
#' @export
#' @inheritSection RAMClustR::ramclustR @return
#' @usage
#' clustering(
#' xcmsObj = NULL,
#' ms = NULL,
#' idmsms = NULL,
#' taglocation = "filepaths",
#' MStag = NULL,
#' idMSMStag = NULL,
#' featdelim = "_",
#' timepos = 2,
#' st = NULL,
#' sr = NULL,
#' maxt = NULL,
#' deepSplit = FALSE,
#' blocksize = 2000,
#' mult = 5,
#' hmax = NULL,
#' sampNameCol = 1,
#' collapse = TRUE,
#' usePheno = TRUE,
#' mspout = TRUE,
#' ExpDes = NULL,
#' normalize = "TIC",
#' qc.inj.range = 20,
#' order = NULL,
#' batch = NULL,
#' qc = NULL,
#' minModuleSize = 2,
#' linkage = "average",
#' mzdec = 3,
#' cor.method = "pearson",
#' rt.only.low.n = TRUE,
#' fftempdir = NULL,
#' replace.zeros = TRUE
#' )
#'
clustering <- function(...){
  RC <- RAMClustR::ramclustR(...)
  RC
}

#' Cluster annotation
#'
#' Cluster annotation using the `InterpretMSSpectrum::findMain`. It creates a
#' pseudo-spectrum with the features grouped in each cluster within a given mz
#' error. The algorithm proposes a representative ion of each
#' cluster/pseudo-spectrum for data reduction.
#'
#' @inheritParams RAMClustR::do.findmain
#' @inherit RAMClustR::do.findmain
#'
#' @inheritSection RAMClustR::do.findmain @return
#' @usage
#' do.findmain(
#' ramclustObj = NULL,
#' cmpd = NULL,
#' mode = "positive",
#' mzabs.error = 0.005,
#' ppm.error = 10,
#' ads = NULL,
#' nls = NULL,
#' scoring = "auto",
#' plot.findmain = TRUE,
#' writeMat = TRUE,
#' writeMS = TRUE,
#' use.z = TRUE
#' )
#' @export
#'
do.findmain <- function(...){
  RC <- RAMClustR::do.findmain(...)
  RC
}


#' Labelled adducts
#'
#' Function ton extract labelled adducts from [clustering].
#'
#' @param RC hclust object from the [do.findmain] function.
#'
#' @return A list containing data frames with (1) all labelled adducts, (2)
#'   representative ion in each cluster/pseudospectrum, and (3) the clustering
#'   data.
#' @export
#'
labelled_adducts <- function (RC) {
# Selection of max intensity ion as cluster representative
Max_int<- lapply(RC$M.ann, function(x) x[which.max(x$int), ])
representative_ions <- dplyr::bind_rows(Max_int)
representative_ions$name <- paste(round(representative_ions$mz,4),
                                  round(RC$clrt,2),
                                  sep = "_")

# Abducts of representative ions
Representative_adducts <- sapply(RC$M.ann, function(x) x[which.max(x$int), ]$adduct)

# Selection of labelled max intensity ion as cluster representative
Labelled_int <- lapply(RC$M.ann, function(x) {
  xl <- x[which(!is.na(x$label)), ]
  xl[which.max(xl$int), ]
})
Labelled_ions <- dplyr::bind_rows(Labelled_int)

# We save most important cluster data
cluster_data <- data.frame(cluster = RC$ann,
                           Max_int_ion_mz = representative_ions$mz,
                           Max_int_ion_adduct = Representative_adducts,
                           Labelled_ion_mz = Labelled_ions$mz,
                           Labelled_ion_adduct = Labelled_ions$label,
                           RC_mz = RC$M,
                           retention_time = RC$clrt)

# All labelled ions
All_labelled_adducts <- lapply(RC$M.ann, function(x) {
  xl <- x[which(!is.na(x$label)), ]
})
All_labelled_adducts <- dplyr::bind_rows(All_labelled_adducts)
return(list(All_labelled_adducts = All_labelled_adducts, Labelled_ions = Labelled_ions, representative_ions = representative_ions, cluster_data = cluster_data))
}

#' Feature reduction
#'
#' `feature_reduction` is a function to reduce the feature table using the info
#' from the [clustering]. A representative feature/ion is selected from each
#' cluster discarding the rest of the feature of the same cluster (normally, the
#' most intense one).
#'
#' @param mdata data matrix with samples in rows and features in columns
#' @param representative_ions data frame from the list given by [labelled_adducts] function
#' @param RC hclust object from the [do.findmain] function
#'
#' @return A reduced data frame with important features based on the clustering
#'   and features no clustered.
#' @export
#' @usage
#' feature_reduction(mdata,
#' representative_ions,
#' RC)
#'
feature_reduction <- function (mdata, representative_ions, RC) {

  mdata <- as.matrix(mdata)

mz <- colnames(mdata) %>%
    stringr::str_split(.,"\\_") %>%
    lapply(.,function(x) x[1]) %>%
    unlist() %>%
    as.numeric()

# Clustered features
xdata_cluster_ions <- mdata[, mz %in% representative_ions$mz]

# Singletons
clustered_mz<- lapply(RC$M.ann,
                      function (x) x$mz) %>% unlist() %>%  as.numeric()

message("A number of ",
        length(clustered_mz),
        " features have been clustered into ",
        dim(xdata_cluster_ions)[2],
        " representative features")

mdata <- as.data.frame(mdata)
xdata_cluster_ions <- as.data.frame(xdata_cluster_ions)

singletons <- mdata[,!mz %in% clustered_mz]

message("A number of ",
        RC$nsing,
        " features correspond to singletons")
#Combine singletons and molecular ions
xdata_reduced <- cbind.data.frame(xdata_cluster_ions, singletons)

message("Original dataset has ",
        ncol(mdata),
        " features")

message("Cluster representative ions dataset has ", ncol(xdata_cluster_ions),
        " features")

message("Singletons dataset has ", ncol(singletons), " features")

message("Reduced dataset has ", ncol(xdata_reduced), " features")
return(xdata_reduced)
}
