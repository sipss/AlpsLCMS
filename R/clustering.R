#' Get st parameter value for RamclustR
#'
#' @param XObj A lcms_dataset
#' @return recomended st parameter value for RamclustR
#' @export
#'
getRamSt <- function(XObj) {
  featInfo <- featureDefinitions(XObj)
  histo <- hist((featInfo$rtmax-featInfo$rtmin)/2)
  st <- round(median(featInfo$rtmax-featInfo$rtmin)/2,
              digits = 2)
  abline(v=st)
  return(st)
}

#' ramclustR
#'
#' Main clustering function for grouping features based on their analytical behavior.
#'
#' @export
#'
ramclustR <- function(ms,
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

#' do.findmain
#'
#' Cluster annotation function: inference of 'M' - molecular weight of the compound giving rise to each spectrum - using the InterpretMSSpectrum::findMain function
#'
#' @export
#'
do.findmain <- function(RC,
                   nls,
                   mode,
                   mzabs.error,
                   ppm.error,
                   writeMat){
  RC_fm <- RAMClustR::do.findmain <- function(RC = RC,
                          nls = nls,
                          mode = mode,
                          mzabs.error = mzabs.error,
                          ppm.error = ppm.error,
                          writeMat = writeMat)
    return(RC_fm)
}
