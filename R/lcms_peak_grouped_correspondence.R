#' Peak Correspondence
#'
#' Peak correspondence is carried out by the *'lcms_groupChromPeaks'* method,
#' with parameters obtained form `IPO`. Peak Correspondece consist in
#' grouping peaks on retention time axis with the purspose of associate
#' them to spectra on the mass/chage axis.
#' #' Note: signal processing generally  consists in three main steps:
#' (1) peak detection (`lcmfindChromPeaks_cwp` function),
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
#' \dontrun{
#' file_name <-  system.file("extdata", "peakdet_align.rds", package = "NIHSlcms")
#' peakdet_align <- base::readRDS(file_name)
#' opt_result_path <-  system.file("extdata", package = "NIHSlcms")
#' lcms_preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#'
#' peak_table <- lcms_group_peaks(peakdet_align, params = lcms_preproc_params)
#' print(peak_table)
#' }
lcms_group_peaks <- function (peakdet_align, params) {
      quiet <- function(x) {
        base::sink(base::tempfile())
        base::on.exit(base::sink())
        base::invisible(base::force(x))
        }

      cat("\n","Grouping peaks using the optimized set of parameters obtained from IPO package.","\n")
  pdp <- base::suppressWarnings(
            base::suppressMessages(quiet(xcms::PeakDensityParam(sampleGroups = peakdet_align$treatment,
                                                                binSize = params$mzwid,
                                                                minFraction = params$minFraction,
                                                                minSamples = params$minSamples,
                                                                maxFeatures =params$maxFeatures,
                                                                bw = params$bw))
                                     )
                              )

  peak_table <- base::suppressWarnings(
                    base::suppressMessages(quiet(xcms::groupChromPeaks(peakdet_align, param = pdp))
                                           )
                                      )
  return(peak_table)
}
