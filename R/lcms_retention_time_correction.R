#' Retention Time Correction
#'
#' Retention time correction is performed using *'obiwarp'* method.
#' Its optimum parameters are obtaine from `IPO` Package.
#' Note: signal processing generally  consists in three main steps:
#' (1) peak detection (`lcms_find_chrom_peaks_cwp` function),
#' (2) peak alignment (`lcms_align_rtime` function) and
#' (3) peak correspondence (`lcms_group_peaks` function).
#'
#' @param peakdet A dataset with detected peaks from the
#' `lcms_find_chrom_peaks_cwp` function
#' @param params A converted parameters template from IPO parameters.
#' @examples
#' \dontrun{
#' file_name <-  system.file("extdata", "peakdet.rds", package = "NIHSlcms")
#' peakdet <- base::readRDS(file_name)
#' opt_result_path <-  system.file("extdata", package = "NIHSlcms")
#' preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#'
#' peakdet_align <- lcms_align_rtime(peakdet, params = preproc_params)
#' print(peakdet_align)
#' }
#'
#' @return A dataset with (1) detected (`lcms_find_chrom_peaks_cwp` function)
#' and (2) aligned (`lcms_align_rtime` function) peaks
#' @export
lcms_align_rtime <- function (peakdet, params) {

  quiet <- function(x) {
      base::sink(base::tempfile())
      base::on.exit(base::sink())
      base::invisible(base::force(x))
    }

  cat("\n","Aligning the peak table using the optimized set of parameters obtained from IPO package.","\n")

  obiwarp_params <- base::suppressWarnings(
                            base::suppressMessages(
                                    quiet(xcms::ObiwarpParam(binSize = params$profStep,
                                     centerSample = params$centerSample,
                                     response = params$response,
                                     distFun = params$distFun,
                                     gapInit = params$gapInit,
                                     gapExtend = params$gapExtend,
                                     factorDiag = params$factorDiag,
                                     factorGap = params$factorGap,
                                     localAlignment = params$localAlignment,
                                     initPenalty = params$initPenalty))
                                                  )
                                           )

  peakdet_align <-  base::suppressWarnings(
                            base::suppressMessages(
                                    quiet(xcms::adjustRtime(peakdet, param = obiwarp_params))
                                                  )
                                          )
  return(peakdet_align)
}
