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
#' Cluster annotation using the [InterpretMSSpectrum::findMain].
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


featureReduction <- function (RC) {
# Selection of max intensity ion as cluster representative
Max_int<- lapply(RC$M.ann, function(x) x[which.max(x$int), ])
Representative_ions <- dplyr::bind_rows(Max_int)
Representative_ions$name <- paste(round(Representative_ions$mz,4),
                                  round(RC$clrt,2),
                                  sep = "_")

# Abducts of representative ions
Representative_adducts <- sapply(RC$M.ann, function(x) x[which.max(x$int), ]$adduct)

# Selection of labeled max intensity ion as cluster representative
Labeled_int <- lapply(RC$M.ann, function(x) {
  xl <- x[which(!is.na(x$label)), ]
  xl[which.max(xl$int), ]
})
Labeled_ions <- dplyr::bind_rows(Labeled_int)

# We save most important cluster data
cluster_data <- data.frame(cluster = RC$ann,
                           Max_int_ion_mz = Representative_ions$mz,
                           Max_int_ion_adduct = Representative_adducts,
                           labeled_ion_mz = Labeled_ions$mz,
                           labeled_ion_adduct = Labeled_ions$label,
                           RC_mz = RC$M,
                           retention_time = RC$clrt)

# All labeled ions
All_labeled_adducts <- lapply(RC$M.ann, function(x) {
  xl <- x[which(!is.na(x$label)), ]
})
All_labeled_adducts <- dplyr::bind_rows(All_labeled_adducts)
}
