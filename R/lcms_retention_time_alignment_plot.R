#' Retention Time Correction Plot
#'
#' It plots the retention time correction vs the original retention time for each of the samples
#' coloured by sample class.
#'
#' @param data An aligned lcms_dataset
#' @examples
#' file_name <-  system.file("extdata", "peakdet_align.rds", package = "NIHSlcms")
#' data <- base::readRDS(file_name)
#' rta_plot <- lcms_retention_time_alignment_plot(data)
#' print(rta_plot)
#' @return The plot for the retention time correction
#' @export
#'
lcms_retention_time_alignment_plot <- function (data){
  min2sec <- 60
  rtc_df <-  tibble::tibble(
    fileIdx = Biobase::fData(data)$fileIdx,
    treatment = Biobase::pData(data)$treatment[fileIdx],
    ret_time_adj = MSnbase::rtime(data, adjust = TRUE),
    ret_time_orig = MSnbase::rtime(data, adjust = FALSE)
  )
  rt <- round(base::range(rtc_df$ret_time_orig) / min2sec)
  tick_values <- seq(from = rt[1], to = rt[2], by = 2)

  rta_plot <- ggplot2::ggplot(rtc_df) +
    ggplot2::geom_line( ggplot2::aes(x = ret_time_orig / min2sec, y = (ret_time_adj-ret_time_orig), color = treatment, group = fileIdx)) +
    ggplot2::scale_x_continuous("Original retention time (min.)", limits = rt, breaks = tick_values) +
    ggplot2::scale_y_continuous("Retention time correction (s)") +
    ggplot2::ggtitle("Retention time alignment warping for each sample")
}
