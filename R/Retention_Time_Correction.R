#' Retention Time Correction
#'
#' Retention time correction is performed  using *'obiwarp'* method.
#' Its optimum parameters are obtaine from `IPO` Package.
#'
#' @param peakdet
#' @param params A converted parameters template from IPO parameters.
#' @example
#' params <- convert_IPO_to_XCMS(IPO_params)
#' peakdet <- findChromPeaks_cwp(dataset, params = params)
#' peakdet_align <- align_Rtime(peakdet, params = params)
#'
#' @return
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
