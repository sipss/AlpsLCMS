#' Principal Component Analysis (PCA)
#'
#' It performs PCA using on the annotated peak table obtained from a MAIT object.
#'
#' @param MAIT.object MAIT object where it is found an annotated peak table
#' @param treament_col Treatment for the samples.
#' @param Log Set to TRUE if the data should be plotted using the logarithm of the intensity.
#' @param Set to TRUE if the data should be centered around its mean. See scale.
#' @param scale Set to TRUE if the data should be scaled.

#' @return A MAIT.object. Additionaly: Three different PCA scoreplots are printed in three png files.
#' One using PC1 vs PC2, another with PC1 vs PC3 and the last one with PC2 vs PC3.
#' The files will be stored in the directory /PCA_Scoreplots.
#' @export
#' @examples
#'\dontrun{
#' file_name_1 <-  system.file("extdata","peak_table_sig_ANN.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name_1)
#' file_name_2 <-  system.file("extdata","lcms_dataset_pos_rt_rs.rds", package = "NIHSlcms")
#' lcms_dataset <-  base::readRDS(file_name_2)
#' treatment_col <- scales::hue_pal()(length(unique(lcms_dataset$treatment)))

#' peak_table_pca <- lcms_peak_table_pca(peak_table,
#'                                       treament_col = treatment_col,
#'                                       Log = FALSE,
#'                                       center = TRUE, scale = FALSE)
#' print(peak_table_pca)
#' }


lcms_peak_table_pca <- function (MAIT.object = NULL,treament_col, Log = FALSE, center = TRUE, scale = TRUE)
{
  if (is.null(MAIT.object)) {
    stop("No input MAIT object file was given")
  }
    if (is.null(treament_col)) {
    stop("No input treatment column was given")
  }
  if (length(featureSigID(MAIT.object)) == 0) {
    stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and peakAggregation were launched")
  }
  parameters <- list(Log, center, scale)
  names(parameters) <- c("PCA data logarithm", "PCA data centered",
    "PCA data scaled")
  MAIT.object@RawData@parameters@plotPCA <- parameters
  lcms_writeParameterTable(parameters(MAIT.object), folder = MAIT.object@PhenoData@resultsPath)
  data <- scores(MAIT.object)
  clases <- classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  xsaFA <- MAIT.object@RawData@data
  resultsPath <- MAIT.object@PhenoData@resultsPath
  index <- featureSigID(MAIT.object)
  cols <- matrix(nrow = 1)
  textCols <- matrix(nrow = 1)
  for (i in 1:length(clases)) {
    cols <- c(cols, rep(i, classNum[i]))
  }
  cols <- as.character(cols[-1])
  textCols <- 1:length(clases)
  if (Log == FALSE) {
    data <- (scale(t(data[index, ]), center = center, scale = scale))
  }else {
    data <- (scale(t(log10(data[index, ] + 1)), center = center,
      scale = scale))
  }
  if (!file.exists(paste(resultsPath, "PCA", sep = "/"))) {
    dir.create(paste(resultsPath, "PCA", sep = "/"))
  }else {
    cat(" ", fill = TRUE)
    cat(paste("Warning: Folder", paste(resultsPath, "PCA_Results",
      sep = "/"), "already exists. Possible file overwritting.",
      sep = " "), fill = TRUE)
  }


  model <- prcomp(data)
  var_expl <-100 *((model$sdev ^ 2)/sum(model$sdev ^ 2))
  model_df <- data.frame(pc = seq(1, dim(model$x)[1], by = 1), model$x,
                         var_expl = var_expl,
                         treatment = clases,
                         stringsAsFactors = FALSE)


  var_plot <- ggplot2::ggplot(model_df) +
              ggplot2::geom_point(ggplot2::aes(x = pc , y = var_expl), color = "blue") +
              ggplot2::scale_x_continuous("Principal Component") +
              ggplot2::scale_y_continuous("Percentage of Variance (%)") +
              ggplot2::ggtitle(paste("Percentage of Variance per Principal Component"))
  print(var_plot)
  base::suppressWarnings(
    base::suppressMessages(ggplot2::ggsave(paste(paste(resultsPath,
                                              "PCA/Percentage_of_Variance_per_PC.png",
                                              sep = "/")))))

  PC1_PC2_plot  <- ggplot2::ggplot(model_df) +
                   ggplot2::geom_point(ggplot2::aes(x = PC1, y = PC2, color = treatment)) +
                   ggplot2::geom_hline(yintercept = 0, linetype="dashed", color = "black") +
                   ggplot2::geom_vline(xintercept = 0, linetype="dashed", color = "black") +
                   ggplot2::scale_color_manual("Treatment", values = treatment_col) +
                   ggplot2::scale_x_continuous(paste0("PC1 (", round(model_df$var_expl[1], 2), " %)")) +
                   ggplot2::scale_y_continuous(paste0("PC2 (", round(model_df$var_expl[2], 2), " %)"))  +
                   ggplot2::ggtitle(paste("PCA Scoreplot", "PC1", "vs", "PC2"))
  print(PC1_PC2_plot)
   base::suppressWarnings(
    base::suppressMessages(ggplot2::ggsave(paste(paste(resultsPath,"PCA/Scoreplot_PC12.png", sep = "/")),
         plot = PC1_PC2_plot)))




  PC1_PC3_plot  <- ggplot2::ggplot(model_df) +
  ggplot2::geom_point(ggplot2::aes(x = PC1, y = PC3, color = treatment)) +
  ggplot2::geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  ggplot2::geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  ggplot2::scale_color_manual("Treatment", values = treatment_col) +
  ggplot2::scale_x_continuous(paste0("PC1 (", round(model_df$var_expl[1], 2), " %)")) +
  ggplot2::scale_y_continuous(paste0("PC3 (", round(model_df$var_expl[3], 2), " %)"))  +
  ggplot2::ggtitle(paste("PCA Scoreplot", "PC1", "vs", "PC3"))
  print(PC1_PC3_plot)
  base::suppressWarnings(
    base::suppressMessages(ggplot2::ggsave(paste(paste(resultsPath,
                                              "PCA/Scoreplot_PC13.png", sep = "/")))))


  PC2_PC3_plot  <-ggplot2::ggplot(model_df) +
  ggplot2::geom_point(ggplot2::aes(x = PC2, y = PC3, color = treatment)) +
  ggplot2::geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  ggplot2::geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  ggplot2::scale_color_manual("Treatment", values = treatment_col) +
  ggplot2::scale_x_continuous(paste0("PC1 (", round(model_df$var_expl[2], 2), " %)")) +
  ggplot2::scale_y_continuous(paste0("PC2 (", round(model_df$var_expl[3], 2), " %)"))  +
  ggplot2::ggtitle(paste("PCA Scoreplot", "PC2", "vs", "PC3"))
  print(PC2_PC3_plot)
  base::suppressWarnings(
    base::suppressMessages(ggplot2::ggsave(paste(paste(resultsPath,
                                              "PCA/Scoreplot_PC23.png", sep = "/")))))

 loadings_df <- data.frame(model$rotation[, 1:3],
                           feature = 1:dim(model$rotation[, 1:3])[1]) %>%
                           tidyr::gather(key = "PC", value = "loadings", -feature)


 loadings_plot  <- ggplot2::ggplot(loadings_df) +
                   ggplot2::geom_line(ggplot2::aes(x = feature, y = loadings, color = PC)) +
                   ggplot2::ggtitle(paste("Loadings of the PCA model (3 PCs)")) +
                   ggplot2::scale_x_continuous("Signicant Feature Index") +
                   ggplot2::scale_y_continuous("Loadings (adim.)")
 base::suppressWarnings(
   base::suppressMessages(ggplot2::ggsave(paste(paste(resultsPath,
                                                      "PCA/Loadings.png", sep = "/")))))
 print(loadings_plot)

  MAIT.object@FeatureData@pcaModel <- list(model)
  return(MAIT.object)
}

