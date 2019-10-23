#' Base peak chromatogram
#'
#' Base peak chromatograms with retention time axis in minutes.
#'
#' @param chromatogram  A XChromatograms object
#' @param treatment_col Color code by groups. They can be generated using
#' treatment_col <- scales::hue_pal()(length(unique(lcms_dataset$treatment)))
#' names(treatment_col) <- unique(lcms_dataset$treatment)
#' @param rtlim retention time boundaries (e.g. c(4,14))
#'
#' @return A base peak chromatogram
#' @export
#'
#' @examples
#'
#' file_name <- system.file("extdata",
#'                          "lcms_dataset_pos.rds",
#'                           package = "NIHSlcms")
#' lcms_dataset <- lcms_dataset_load(file_name)
#' treatment_col <- scales::hue_pal()(length(unique(lcms_dataset$treatment)))
#' names(treatment_col) <- unique(lcms_dataset$treatment)
#' base_peaks <- xcms::chromatogram(lcms_dataset, aggregationFun = "max")
#' lcms_plotChrom(base_peaks, treatment_col, rtlim = c(4, 14))
#'
#'
lcms_plotChrom <- function (chromatogram, treatment_col, rtlim = NULL){
  min2sec <- 60
  message("Make sure that the column that contains the group class is called `treatment`")

  #we need to modify this in order to be more flexible (treatment_col)
  ret_times <- lapply(chromatogram, FUN = rtime)
  intensities <- lapply(chromatogram, FUN = intensity)

  plot(ret_times[[1]] / min2sec, intensities[[1]], type = "l",
       col = treatment_col[chromatogram$treatment][1], lwd = 1,
       xlab = "Retention time (min)", ylab = "Intensity  (A.U.)",
       xlim = rtlim,
       main = "Base Peak Chromatogram")

  for (i in 2:length(ret_times)){ # we need to modify the sequence using seq_along
    points(ret_times[[i]] / min2sec, intensities[[i]], lwd = 1,
           xlim = rtlim,
           type= "l", col = treatment_col[chromatogram$treatment][i])
    legend("topright", legend = names(treatment_col), fill = treatment_col)
  }
}
