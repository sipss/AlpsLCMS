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
  graphics::abline(v=st)
  return(st)
}

#' clustering of peaks
#'
#' Main clustering function for grouping features based on their analytical
#' behavior.
#' @inheritParams RAMClustR::ramclustR
#' @export
#'
lcms_clustering <- function(ms,
          featdelim,
          st,
          sr,
          ExpDes,
          normalize,
          sampNameCol,
          fftempdir){
  RC <- RAMClustR::ramclustR(ms = ms,
            featdelim = featdelim,
            st = st,
            sr = sr,
            ExpDes = ExpDes,
            normalize = normalize,
            sampNameCol = sampNameCol,
            fftempdir = fftempdir)
  return(RC)
}

#' do_findmain
#'
#' Cluster annotation function: inference of 'M' - molecular weight of the
#' compound giving rise to each spectrum - using the
#' InterpretMSSpectrum::findMain function
#'
#' @export
#' @inheritParams RAMClustR::do.findmain
do_findmain <- function(ramclustObj = NULL,
                        cmpd = NULL,
                        mode = "positive",
                        mzabs.error = 0.005,
                        ppm.error = 10,
                        ads = NULL,
                        nls = NULL,
                        scoring = "auto",
                        plot.findmain = TRUE,
                        writeMat = TRUE,
                        writeMS = TRUE,
                        use.z = TRUE){
  RC_fm <- RAMClustR::do.findmain (ramclustObj = ramclustObj,
                                   cmpd = cmpd,
                                   mode = mode,
                                   mzabs.error = mzabs.error,
                                   ppm.error = ppm.error,
                                   ads = ads,
                                   nls = nls,
                                   scoring = scoring,
                                   plot.findmain = plot.findmain,
                                   writeMat = writeMat,
                                   writeMS = writeMS,
                                   use.z = use.z)
    return(RC_fm)
}
