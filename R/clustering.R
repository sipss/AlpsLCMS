#' Get st parameter value for RamclustR
#'
#' @param XObj A lcms_dataset
#' @return recomended st parameter value for RamclustR
#' @export
#'
getRamSt <- function(XObj) {
  featInfo <- xcms::featureDefinitions(XObj)
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
#' retention time similarity `st` and for correlational similarity.
#'
#' @inheritDotParams RAMClustR::ramclustR
#' @inherit RAMClustR::ramclustR
#'
#' @return hclust object
#' @export
#'
clustering <- function(...){
  RC <- RAMClustR::ramclustR(...)
  RC
}

#' Cluster annotation
#'
#' Cluster annotation using the `InterpretMSSpectrum::findMain`.
#'
#' @inheritDotParams RAMClustR::do.findmain
#' @inherit RAMClustR::do.findmain
#'
#' @return hclust object
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
#' @param RC hclust object from the [do.findmain] function
#'
#' @return
#' @export
#'
labelled_adducts <- function (RC) {
# Selection of max intensity ion as cluster representative
Max_int<- lapply(RC$M.ann, function(x) x[which.max(x$int), ])
Representative_ions <- dplyr::bind_rows(Max_int)
Representative_ions$name <- paste(round(Representative_ions$mz,4),
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
                           Max_int_ion_mz = Representative_ions$mz,
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
return(list(All_labelled_adducts = All_labelled_adducts, Labelled_ions = Labelled_ions, Representative_ions = Representative_ions, cluster_data = cluster_data))
}


#' Feature reduction
#'
#' Function to perform feature reduction from the [clustering]
#'
#' @param mdataImputed data matrix with samples in rows and features in columns
#' @param Representative_ions dataframe from the list given by [labelled_adducts] function
#' @param RC hclust object from the [do.findmain] function
#'
#' @return
#' @export
#'
feature_reduction <- function (mdataImputed, Representative_ions, RC) {

  mdataImputed <- as.matrix(mdataImputed)

# Clustered features
xdata_cluster_ions <- mdataImputed[, mz %in% Representative_ions$mz]

# Singletons
clustered_mz<- lapply(RC$M.ann,
                      function (x) x$mz) %>% unlist() %>%  as.numeric()

message("A number of ",
        length(clustered_mz),
        " features have been clustered into ",
        dim(xdata_cluster_ions)[2],
        " representative features")

mdataImputed <- as.data.frame(mdataImputed)
xdata_cluster_ions <- as.data.frame(xdata_cluster_ions)

singletons <- mdataImputed[,!mz %in% clustered_mz]

message("A number of ",
        RC$nsing,
        " features correspond to singletons")
#Combine singletons and molecular ions
xdata_reduced <- cbind.data.frame(xdata_cluster_ions, singletons)

message("Original dataset has ",
        ncol(mdataImputed),
        " features")

message("Cluster representative ions dataset has ", ncol(xdata_cluster_ions),
        " features")

message("Singletons dataset has ", ncol(singletons), " features")

message("Reduced dataset has ", ncol(xdata_reduced), " features")
return(xdata_reduced)
}
