#' Total Ion Count (TIC)
#'
#' The function performs the Total Ion Count (TIC) for the polarity samples.
#' Function `lcms_tics` stores summarizes TIC information.
#' NOTE: `lcms_tics` assumes that data is already filtered by polarity.
#' @param lcms_dataset An [lcms_dataset_family] object
#' @param treatment Class groups of the samples
#' @return Total Ion Count (TIC) for the polarity samples.
#' @export
#' @family lcms_dataset functions
#' @family lcms_dataset_peak_table functions
lcms_tics <- function(lcms_dataset, treatment = treatment){
  tics <- tibble::tibble(
    file = MSnbase::fromFile(lcms_dataset),
    fileName = Biobase::pData(lcms_dataset)$sampleNames[file],
    treatment = Biobase::pData(lcms_dataset)$treatment[file],
    ret_time = MSnbase::rtime(lcms_dataset),
    polarity = rep(unique(MSnbase::polarity(lcms_dataset),length(file))),
    tic = MSnbase::tic(lcms_dataset)
  )
  #Files sorted by treatment
  tics$fileName <- factor(tics$fileName,
                          levels = unique(tics$fileName)
                          [base::order(lcms_dataset$treatment)])
  return(tics)
}


#' Total Ion Count (TIC) plot
#'
#' The function performs a plot with the Total Ion Count (TIC).
#' @param tics A Total Ion Count object generated with `lcms_tics`
#' @param treatment Class groups of the samples
#' @param rt Retention time
#' @param plot_type The plot class, either boxplot or spectra
#' @return Total Ion Count (TIC) for the polarity samples.
#' @export
#' @family lcms_dataset functions
#' @family lcms_dataset_peak_table functions
lcms_plot_tics <- function(tics, treatment = treatment, rt = NULL, plot_type = "spec"){
  min2sec <- 60
  treatment_col <- scales::hue_pal()(length(unique(tics$treatment)))
  names(treatment_col) <- unique(tics$treatment)

  if (unique(tics$polarity) == 1){
    polarity <- "(positive polarity)"
  } else if (unique(tics$polarity) == 0){
    polarity <- "(negative polarity)"
  }

  if(is.null(rt)){
    rt <-  round (base::range(tics$ret_time) /  min2sec)
  } else{
    rt <- round(rt)
  }

  if(plot_type == "spec"){
    if(diff(rt)> 6){
      tick_values <- seq(from = rt[1], to = rt[2], by = 2)
    }else if (diff(rt) <= 6 & diff(rt) >= 1){
      tick_values <- seq(from = rt[1], to = rt[2], by = 1)
    }else {
      stop("Non valid retention time range")
    }
    ggplot2::ggplot(tics) +
      ggplot2::geom_line(ggplot2::aes(x = ret_time /  min2sec, y = tic, color = treatment, group = file)) +
      ggplot2::scale_x_continuous("Retention time (min)", limits = rt, breaks = tick_values) +
      ggplot2::scale_y_continuous("Total Ion Count (a.u.)") +
      ggplot2::scale_colour_manual("Treatment", values = treatment_col) +
      ggplot2::ggtitle(paste("Total Ion Count across all retention times", polarity))
  }else if(plot_type == "boxplot"){
    tics <- tics %>% dplyr::filter(ret_time >= rt[1] * min2sec & ret_time <= rt[2]* min2sec)
    ggplot2::ggplot(tics) +
      ggplot2::geom_boxplot(ggplot2::aes(x = fileName, y = tic, fill = treatment)) +
      ggplot2::scale_fill_manual("Treatment", values = treatment_col) +
      ggplot2::scale_y_log10() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::ggtitle(paste("Boxplot of the Total Ion Count by sample", polarity))
  }
}
