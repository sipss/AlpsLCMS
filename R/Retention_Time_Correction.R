#' Retention Time Correction
#'
#' Retention time correction is performed using *'obiwarp'* method.
#' Its optimum parameters are obtaine from `IPO` Package.
#' Note: signal processing generally  consists in three main steps:
#' (1) peak detection (`findChromPeaks_cwp` function),
#' (2) peak alignment (`align_Rtime` function) and
#' (3) peak correspondence (`group_peaks` function).
#'
#' @param peakdet A lcms_dataset with detected peaks from the
#' `findChromPeaks_cwp` function
#' @param params A converted parameters template from IPO parameters.
#' @examples
#' params <- convert_IPO_to_XCMS(IPO_params)
#' peakdet <- findChromPeaks_cwp(dataset, params = params)
#' peakdet_align <- align_Rtime(peakdet, params = params)
#'
#' @return A lcms_dataset with (1) detected (`findChromPeaks_cwp` function)
#' and (2) aligned (`align_Rtime` function) peaks
#' @export
#'
align_Rtime <- function (peakdet, params) {
  obiwarp_params <- xcms::ObiwarpParam(binSize = params$profStep,
                                     centerSample = params$centerSample,
                                     response = params$response,
                                     distFun = params$distFun,
                                     gapInit = params$gapInit,
                                     gapExtend = params$gapExtend,
                                     factorDiag = params$factorDiag,
                                     factorGap = params$factorGap,
                                     localAlignment = params$localAlignment,
                                     initPenalty = params$initPenalty)

  peakdet_align <- xcms::adjustRtime(peakdet, param = obiwarp_params)
  return(peakdet_align)
}
