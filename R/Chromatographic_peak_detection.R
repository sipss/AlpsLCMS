#' Chromatographic peak detection (CentWave)
#'
#' The findChromPeaks_cwp function performs the chromatographic peak
#' detection on LC/GC-MS data. The standard method for peak detection
#' is *'CentWave'*. We must initialize its parameters according to the
#' `IPO` Package optimization. Peak detection aims to detect important
#' features (peaks) on the chromatographic axis. This will be useful
#' for a posterior peak alignment on the chormatophic axis.
#'
#' @param lcms_dataset An [lcms_dataset_family] object
#' @param params A converted parameters template from IPO parameters.
#' @example
#' params <- convert_IPO_to_XCMS(IPO_params)
#' peakdet <- findChromPeaks_cwp (dataset, params = params)
#'
#' @return
#' @export
#'
findChromPeaks_cwp <- function (lcms_dataset, params) {
  cwp <- xcms::CentWaveParam(peakwidth = params$peakwidth,
                           ppm = params$ppm,
                           mzdiff = params$mzdiff,
                           snthresh =  params$snthresh,
                           noise = params$noise,
                           prefilter = params$prefilter,
                           mzCenterFun = params$mzCenterFun,
                           integrate = params$integrate,
                           fitgauss = params$fitgauss,
                           verboseColumns = params$verbose.columns)

  peakdet <- xcms::findChromPeaks(dataset, param = cwp)
  return(peakdet)
}


