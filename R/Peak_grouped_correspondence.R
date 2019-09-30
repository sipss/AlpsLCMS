#' Peak Correspondence
#'
#' Peak correspondence is carried out by the *'groupChromPeaks'* method,
#' with parameters obtained form `IPO`. Peak Correspondece consist in
#' grouping peaks on retention time axis with the purspose of associate
#' them to spectra on the mass/chage axis.
#' #' Note: signal processing generally  consists in three main steps:
#' (1) peak detection (`findChromPeaks_cwp` function),
#' (2) peak alignment (`align_Rtime` function) and
#' (3) peak correspondence (`group_peaks` function).
#' After this stage the peak table is finally obtained.
#'
#' @param peakdet_align A lcms_dataset with (1) detected (`findChromPeaks_cwp`
#' function) and (2) aligned (`align_Rtime` function) peaks
#' @param params A converted parameters template from IPO parameters.
#'
#' @return A lcms_dataset with (1) detected (`findChromPeaks_cwp` function),
#' (2) aligned (`align_Rtime` function) and (3) grouped  (`group_peaks`
#' function) peaks.
#' @export
#'
#' @examples
#' params <- convert_IPO_to_XCMS(IPO_params)
#' peakdet <- findChromPeaks_cwp(dataset, params = params)
#' peakdet_align <- align_Rtime(peakdet, params = params)
#' peak_table <- group_peaks(peakdet_align, params = params)
#'
group_peaks <- function (peakdet_align, params) {
  pdp <- xcms::PeakDensityParam(sampleGroups = peakdet_align$treatment,
                              binSize = params$mzwid,
                              minFraction = params$minFraction,
                              minSamples = params$minSamples,
                              maxFeatures =params$maxFeatures,
                              bw = params$bw)

  peak_table <- xcms::groupChromPeaks(peakdet_align, param = pdp)
  return(peak_table)
}
