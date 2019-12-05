#' Base peak chromatogram
#'
#' Base peak chromatograms with retention time axis in minutes.
#'
#' @param chromatogram_object  A XChromatograms object
#' @param treatment_col Color code by groups.
#' @param rtlim retention time boundaries (e.g. c(4,8))
#' @return A base peak chromatogram
#' @export
#' @examples
#' file_name <- system.file("extdata",
#'                          "dataset_pos.rds",
#'                           package = "NIHSlcms")
#' dataset <- lcms_dataset_load(file_name)
#' treatment_col <- scales::hue_pal()(length(unique(dataset$treatment)))
#' names(treatment_col) <- unique(dataset$treatment)
#' base_peaks <- xcms::chromatogram(dataset, aggregationFun = "max")
#' lcms_plot_chrom(base_peaks, treatment_col, rtlim = c(4, 8))
#'
lcms_plot_chrom <- function (chromatogram_object, treatment_col, rtlim = NULL){
  min2sec <- 60
  message("Make sure that the column that contains the group class is called `treatment`")

  #we need to modify this in order to be more flexible (treatment_col)
  ret_times <- lapply(chromatogram_object, FUN = rtime)
  intensities <- lapply(chromatogram_object, FUN = intensity)

  plot(ret_times[[1]] / min2sec, intensities[[1]], type = "l",
       col = treatment_col[chromatogram_object$treatment][1], lwd = 1,
       xlab = "Retention time (min)", ylab = "Intensity  (A.U.)",
       xlim = rtlim,
       main = "Base Peak Chromatogram")

  for (i in 2:length(ret_times)){ # we need to modify the sequence using seq_along
    points(ret_times[[i]] / min2sec, intensities[[i]], lwd = 1,
           xlim = rtlim,
           type= "l", col = treatment_col[chromatogram_object$treatment][i])
    legend("topright", legend = names(treatment_col), fill = treatment_col)
  }
}
