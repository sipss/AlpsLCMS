#' Chromatographic peak detection (CentWave)
#'
#' The findChromPeaks_cwp function performs the chromatographic peak
#' detection on LC/GC-MS data. The standard method for peak detection
#' is *'CentWave'*. We must initialize its parameters according to the
#' `IPO` Package optimization. Peak detection aims to detect important
#' features (peaks) on the chromatographic axis. This will be useful
#' for a posterior peak alignment on the chormatophic axis.
#' Note: signal processing generally  consists in three main steps:
#' (1) peak detection (`findChromPeaks_cwp` function),
#' (2) peak alignment (`align_Rtime` function) and
#' (3) peak correspondence (`group_peaks` function).
#'
#' @param lcms_dataset An [lcms_dataset_family] object
#' @param params A converted parameters template from IPO parameters.
#' @examples
#' \dontrun{
#' params <- convert_IPO_to_XCMS(IPO_params)
#' peakdet <- findChromPeaks_cwp (dataset, params = params)
#' }
#'
#' @return A lcms_dataset with detected peaks
#' @export
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

  peakdet <- xcms::findChromPeaks(lcms_dataset, param = cwp)
  return(peakdet)
}


