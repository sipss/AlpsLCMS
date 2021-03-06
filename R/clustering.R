#' Get st parameter value for clustering
#'
#' `getRamSt` function calculates the optimal value for the parameter `st` the [AlpsLCMS::clustering]
#' function. It uses the `lcms_dataset` after the xcms preprocessing.
#'
#' @param lcms_dataset A lcms_dataset
#' @return the recommended value for the `st` parameter in [AlpsLCMS::clustering]
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
#' @inheritDotParams RAMClustR::ramclustR
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
#' Cluster annotation using the InterpretMSSpectrum::findMain and [RAMClustR::do.findmain]. It creates a
#' pseudo-spectrum with the features grouped in each cluster within a given mz
#' error. The algorithm proposes a representative ion of each
#' cluster/pseudo-spectrum for data reduction.
#'
#' @inherit RAMClustR::do.findmain
#' @inheritDotParams RAMClustR::do.findmain
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
do_findmain <- function(...){
  RC <- RAMClustR::do.findmain(...)
  RC
}

#' Feature reduction
#'
#' `feature_reduction` is a function to reduce the feature table using the info
#' from the [clustering]. A representative feature/ion is selected from each
#' cluster discarding the rest of the feature of the same cluster (normally, the
#' most intense one).
#'
#' @param RC hclust object from the [do_findmain] function
#' @param xdataImp data matrix with samples in rows and features in columns. The
#'   matrix should not include metadata or any columns without LC-MS data. Note
#'   that depending on the origins of the dataset, colnames might have a
#'   character starting with "X" avoiding a number in the first character (see
#'   `loaded = TRUE`)
#' @param loaded if matrix is loaded into R, an "X" will be added to each
#'   colname. In this case, set `loaded = TRUE`.
#' @return A list containing:
#' $singletons: data frame with features that has not been clustered.
#' $xdata_mol_ions: data frame with the representative ions of each cluster
#' $xdata_reduced: data frame with the reduced dataset (singletons +
#' representative ions of each cluster)
#' @export
#' @usage
#' feature_reduction(RC,
#' xdata,
#' loaded = TRUE)
feature_reduction <- function (RC, xdataImp, loaded = FALSE) {
  # Selection of max intensity ion as cluster representative
  Max_int<- lapply(RC$M.ann, function(x) x[which.max(x$int), ])
  Molecular_ions <- dplyr::bind_rows(Max_int)
  Molecular_ions$name <- paste(round(Molecular_ions$mz,4),
                               round(RC$clrt,2),
                               sep = "_")
  # Get mz and rt columns for the feature table
  if (loaded == TRUE){
  mz <- colnames(xdataImp) %>%
    stringr::str_split(.,"\\_") %>%
    lapply(.,function(x) x[1]) %>%
    unlist() %>%
    substring(., 2) %>% # add if xdata was loaded (cols with an X)
    as.numeric()
  } else {
    mz <- colnames(xdataImp) %>%
      stringr::str_split(.,"\\_") %>%
      lapply(.,function(x) x[1]) %>%
      unlist() %>%
      # substring(., 2) %>% # add if xdata was loaded (cols with an X)
      as.numeric()
  }

  mz %in% Molecular_ions$mz %>% table()
  message("length of features: ", length(mz))

  # Clustered features
  xdata_mol_ions <- xdataImp[, mz %in% Molecular_ions$mz]

  # Singletons
  clustered_mz<- lapply(RC$M.ann,
                        function (x) x$mz) %>%
    unlist() %>%
    as.numeric()

  message("A number of ",
          length(clustered_mz),
          " features have been clustered into ",
          dim(xdata_mol_ions)[2],
          " representative features")
  xdataImp <- as.data.frame(xdataImp)
  xdata_mol_ions <- as.data.frame(xdata_mol_ions)
  singletons <- xdataImp[,!mz %in% clustered_mz]#check xdata[,-1]

  message("A number of ",
          ncol(singletons),
          " features correspond to singletons")
  #Combine singletons and molecular ions
  xdata_reduced <- cbind.data.frame(xdata_mol_ions, singletons)

  message("Original dataset has ",
          ncol(xdataImp),
          " features")
  message("")

  message("Molecular ions dataset has ", ncol(xdata_mol_ions),
          " features")
  message("Singletons dataset has ", ncol(singletons), " features")
  message("Reduced dataset has ", ncol(xdata_reduced), " features")

  return(list(singletons = singletons, xdata_mol_ions = xdata_mol_ions, xdata_reduced = xdata_reduced))
}



#   # Abducts of representative ions
#   Representative_adducts <- sapply(RC$M.ann, function(x) x[which.max(x$int), ]$adduct)
#
#   # Selection of labelled max intensity ion as cluster representative
#   Labelled_int <- lapply(RC$M.ann, function(x) {
#     xl <- x[which(!is.na(x$label)), ]
#     xl[which.max(xl$int), ]
#   })
#   labelled_ions <- dplyr::bind_rows(Labelled_int)
#
#   # We save most important cluster data
#   cluster_data <- data.frame(cluster = RC$ann,
#                              Max_int_ion_mz = representative_ions$mz,
#                              Max_int_ion_adduct = Representative_adducts,
#                              Labelled_ion_mz = labelled_ions$mz,
#                              Labelled_ion_adduct = labelled_ions$label,
#                              RC_mz = RC$M,
#                              retention_time = RC$clrt)
#
#   # All labelled ions
#   all_labelled_adducts <- lapply(RC$M.ann, function(x) {
#     xl <- x[which(!is.na(x$label)), ]
#   })
#   all_labelled_adducts <- dplyr::bind_rows(all_labelled_adducts)
#   return(list(all_labelled_adducts = all_labelled_adducts, labelled_ions = labelled_ions, representative_ions = representative_ions, cluster_data = cluster_data))


# feature_reduction <- function (mdata, representative_ions, RC) {
#   mdata <- as.matrix(mdata)
#
#   # Get mz columns for the feature table
#   mz <- colnames(mdata) %>%
#     stringr::str_split(., "\\_") %>%
#     lapply(., function(x)
#       x[1]) %>%
#     unlist() %>%
#     as.numeric()
#
#   # Ramcluster make a round of 4 decimals, we have to round mz too
#   mz <- round(mz, 4)
#   representative_ions_mz <- round(representative_ions$mz, 4)
#
#   # Clustered features
#   xdata_cluster_ions <- mdata[, mz %in% representative_ions_mz]
#
#   # Singletons
#   clustered_mz <- lapply(RC$M.ann,
#                          function (x)
#                            x$mz) %>% unlist() %>%  as.numeric()
#
#   mdata <- as.data.frame(mdata)
#   xdata_cluster_ions <- as.data.frame(xdata_cluster_ions)
#
#   singletons <- mdata[, !mz %in% clustered_mz]
#
#   message("A number of ",
#           RC$nsing,
#           " features correspond to singletons")
#
#   #Combine singletons and molecular ions
#   xdata_reduced <- cbind.data.frame(xdata_cluster_ions, singletons)
#
#   message("Original dataset has ",
#           ncol(mdata),
#           " features")
#
#   message("Cluster representative ions dataset has ",
#           ncol(xdata_cluster_ions),
#           " features")
#
#   message("Singletons dataset has ", ncol(singletons), " features")
#
#   message("Reduced dataset has ", ncol(xdata_reduced), " features")
#   return(xdata_reduced)
# }
