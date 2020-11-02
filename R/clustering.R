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
