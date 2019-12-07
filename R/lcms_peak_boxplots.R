#' Boxplots for significant peak table features
#'
#' It performs boxplots for any of significant features found in a peak table. All the
#' boxplots are stored in a directory /Boxplots/Boxplot_spectra_
#'
#' @param MAIT.object MAIT object where it is found an annotated peak table
#' @param treament_col Treatment for the samples.
#' @return BoxPlots are stored in folders associated to their corresponding spectra. No explicit plot is produced by the device
#' @export
#' @examples
#' file_name_1 <-  system.file("extdata","peak_table_sig_ann.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name_1)
#' file_name_2 <-  system.file("extdata","dataset_pos_rt_rs.rds", package = "NIHSlcms")
#' dataset <-  base::readRDS(file_name_2)
#' treatment_col <- scales::hue_pal()(length(unique(dataset$treatment)))
#' lcms_peak_table_boxplots(peak_table,
#'                          treatment_col = treatment_col)


lcms_peak_table_boxplots <- function (MAIT.object = NULL, treatment_col) {
   if (is.null(treatment_col)) {
    stop("No input treatment column was given")
  }
  if (is.null(MAIT.object)) {
    stop("No input MAIT object file was given")
  }
  if (length(featureSigID(MAIT.object)) == 0) {
    stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and peakAggregation were launched")
  }
  data <- scores(MAIT.object)
  index <- featureSigID(MAIT.object)
  class <- classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath
  clases <- matrix(nrow = 1)
  for (i in c(1:length(class))) {
    clases <- c(clases, rep(class[i], classNum[i]))
  }
  clases <- clases[-1]
  clases <- as.factor(clases)
  aux <- t(data[index, ])
  if (!file.exists(paste(resultsPath, "Boxplots", sep = "/"))) {
    dir.create(paste(resultsPath, "Boxplots", sep = "/"))
  }else {
    cat(" ", fill = TRUE)
    cat(paste("Warning: Folder", paste(resultsPath, "Boxplots",
      sep = "/"), "already exists. Possible file overwritting.",
      sep = " "), fill = TRUE)
  }
  for (i in c(1:length(index))) {
    peaks_df <- data.frame(peaks = aux[,i],
                           treatment = clases,
                           stringsAsFactors = FALSE)
    ggplot2::ggplot(peaks_df) +
    ggplot2::geom_boxplot(ggplot2::aes(x = treatment, y = peaks, fill = treatment)) +
    ggplot2::scale_fill_manual("Treatment", values = treatment_col) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::ggtitle(paste("Boxplot of the Spectra", index[i]))
    base::suppressWarnings(
      base::suppressMessages(ggplot2::ggsave(paste(paste(resultsPath, "Boxplots/Boxplot_spectra_",
      sep = "/"), index[i], ".png", sep = ""))))
  }
}
