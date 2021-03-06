#' Filter an experiment by its polarity
#'
#' This function was merged in the MSnbase package on 2018-12-10. It filter by
#' positive or negative ionization mode.
#'
#' @param lcms_dataset A MSnExp object.
#' @param polarity. The polarity to keep. set `polarity. = 1` for positive
#'   ionization mode. set `polarity. = 0` for negative ionization mode.
#' @family dataset functions
#' @family filtering functions
#' @return A filtered lcms_dataset with the selected polarity.
#' @export
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_read_samples("yourDir/sample1.mzXML", mode = "onDisk")
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' is.positive(dataset_pos)
#' }
lcms_filter_polarity <- function(lcms_dataset, polarity.) {
  if (missing(polarity.)) return(lcms_dataset)
  lcms_dataset <- MSnbase::filterPolarity(lcms_dataset, polarity. = polarity.)
  lcms_dataset
}

#' Filter by mass/charge
#'
#' The function filters the dataset by m/z.
#'
#' @param dataset A lcms_dataset.
#' @param mz  The range of masses to filter.
#' @family dataset functions
#' @family filtering functions
#' @return A filtered lcms_dataset with the selected m/z range.
#' @export
#' @examples
#' \dontrun{
#' mz_range = c(300, 400)
#' dataset <- lcms_dataset_load(system.file("extdata",
#'                                        "dataset_pos.rds",
#'                                        package = "AlpsLCMS"))
#' dataset_pos_mz <-lcms_filter_mz(dataset, mz = mz_range)
#'
#' print(dataset_pos_mz)
#' }

lcms_filter_mz <- function(dataset, mz){
  dataset <- MSnbase::filterMz(dataset, mz)

}


#' Filter by retention time
#'
#' This function converts seconds into minutes to cut and keep
#' a range of the retention time in minutes.
#'
#' @param dataset A lcms_dataset filtered by retention time.
#' @param rt Range of the retention time to keep in minutes.
#' @return A filtered lcms_dataset with the selected retention time range.
#' @family dataset functions
#' @family filtering functions
#' @export
#' @examples
#' \dontrun{
#' rtime_range = c(5,10)
#' dataset <- lcms_dataset_load(system.file("extdata",
#'                                          "dataset_pos.rds",
#'                                          package = "AlpsLCMS"))
#' dataset_pos_rt <-lcms_filter_rt_min(dataset, rt = rtime_range)
#'
#' print(dataset_pos_rt)
#' }
lcms_filter_rt_min <- function (dataset, rt = c(4, 14)){
  min2sec <- 60
  dataset <- MSnbase::filterRt(dataset, rt = rt * min2sec)
}

#' Filter by sample type
#'
#' In a dataset, there are different types of samples for
#' checking purposes. For instance, quality control (QC) samples
#' may be pools of all samples, and blank samples may be composed
#' by the solvent used. Therefore, these samples should be removed
#' before alignment and processing from the main dataset. Use `NULL`
#' if there is not any especial sample in the dataset.
#'
#' The function can distinguish and filter different samples types:
#' * Regular samples (regular_samples).
#' * Blank samples (blanks).
#' * Quality Control samples (QCs).
#' @param dataset A lcms_dataset
#' @param especial_samples A list with the especial samples names.
#' Use `NULL` if there is not any especial sample in the dataset.
#' @return A list containing three lcms_datasets with different type: regular,
#' quality control and blank samples.
#' @family dataset functions
#' @family filtering functions
#' @export
#' @examples
#' \dontrun{
#' dataset <- lcms_dataset_load(system.file("extdata",
#'                                         "dataset_pos_rt.rds",
#'                                         package = "AlpsLCMS"))
#' especial_samples <- list(QC = NULL, blank = NULL)
#' datasets_by_class_type <- lcms_filter_sample_type(dataset, especial_samples)
#' dataset_pos_rt_rs <- datasets_by_class_type$regular_samples
#' dataset_pos_rt_qcs <- datasets_by_class_type$QCs
#' dataset_pos_rt_bks <- datasets_by_class_type$blanks
#'
#' print(dataset_pos_rt_rs)
#' print(dataset_pos_rt_qcs)
#' print(dataset_pos_rt_bks)
#' }

lcms_filter_sample_type <- function(dataset,  especial_samples){

  QC_label <- especial_samples$QC
  blank_label <- especial_samples$blank

  QC_index <- which(dataset$treatment == QC_label)
  if (length(QC_index) == 0){
    QCs <- NULL
  } else{
    QCs <- MSnbase::filterFile(dataset,file = QC_index)
  }


  blank_index <- which(dataset$treatment == blank_label)
  if (length(blank_index) == 0){
    blanks <- NULL
  } else{
    blanks <- MSnbase::filterFile(dataset,file = blank_index)
  }

  sample_index <- which(!(dataset$treatment %in% c(QC_label, blank_label)))
  if (length(sample_index) == 0){
    regular_samples <- NULL
    stop("Your dataset doesn't have any sample not considered Blank or QC sample")
  } else{
    regular_samples <- MSnbase::filterFile(dataset, file = sample_index)
  }

  datasets_by_class_type <- list(regular_samples = regular_samples,
                                 QCs = QCs, blanks = blanks)
}


#' Chromatographic peak detection (CentWave)
#'
#' The findChromPeaks_cwp function performs the chromatographic peak
#' detection on LC/GC-MS data. The standard method for peak detection
#' is *'CentWave'*. Peak detection aims to detect important
#' features (peaks) on the chromatographic axis. This will be useful
#' for a posterior peak alignment on the chormatophic axis.
#'
#' Note: signal processing generally  consists in three main steps: (1) peak
#' detection (`lcms_find_chromp_peaks_cwp` function), (2) peak alignment
#' (`obiAdjust_Rtime` function) and (3) peak correspondence (`IPO_group_peaks`
#' function). The optimized set of parameters for this signal preprocessing can
#' be obatained from `IPO` Package.
#'
#' @param dataset A lcms_dataset.
#' @param params A converted parameters template from IPO parameters.
#' @return A lcms_dataset with the detected peaks added.
#' @family dataset functions
#' @family peak detection functions
#' @export
#' @examples
#' \dontrun{
#' file_path <- system.file("extdata", "rearrange_mait", package = "AlpsLCMS")
#' samples_mzxml <- list.files(file_path, recursive = TRUE, full.names = TRUE)
#' meta_path <- system.file("extdata", "metadata.xlsx", package = "AlpsLCMS")
#' opt_result_path <-  system.file("extdata", package = "AlpsLCMS")
#' preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#'
#' dataset <- suppressWarnings(lcms_read_samples(samples_mzxml, mode = "onDisk"))
#' metadata <- lcms_meta_read(meta_path)
#' dataset_meta <- lcms_meta_add(dataset, metadata, by = "sampleNames")
#'
#' peakdet <- find_peaks_cwp(dataset_meta, params = preproc_params)
#' print(peakdet)
#' }

find_peaks_cwp <- function (dataset, params) {
  quiet <- function(x) {
    base::sink(base::tempfile())
    base::on.exit(base::sink())
    base::invisible(base::force(x))
  }


  cat("\n","Finding chromatographic peaks","\n")
  base::suppressWarnings(
    base::suppressMessages(
      quiet(require(xcms))))
  cwp <- base::suppressWarnings(
    base::suppressMessages(
      quiet(xcms::CentWaveParam(peakwidth = params$peakwidth,
                                ppm = params$ppm,
                                mzdiff = params$mzdiff,
                                snthresh =  params$snthresh,
                                noise = params$noise,
                                prefilter = params$prefilter,
                                mzCenterFun = params$mzCenterFun,
                                integrate = params$integrate,
                                fitgauss = params$fitgauss,
                                verboseColumns = params$verbose.columns))
    )
  )

  peakdet <- base::suppressWarnings(
    base::suppressMessages(
      quiet(xcms::findChromPeaks(dataset, param = cwp))
    )
  )
  return(peakdet)
}

#' Retention Time Correction
#'
#' Retention time correction is performed using *'obiwarp'* method.
#'
#' Note: signal processing generally  consists in three main steps:
#' (1) peak detection (`find_peaks_cwp` function),
#' (2) peak alignment (`obiAdjust_Rtime` function) and
#' (3) peak correspondence (`IPO_group_peaks` function).
#' The optimized set of parameters for this signal preprocessing can be obatained from `IPO` Package.
#'
#' @param peakdet A lcms_dataset with detected peaks obtained from the
#' `find_peaks_cwp` function.
#' @param params A converted parameters template from IPO parameters.
#' @return A lcms_dataset with (1) detected (`find_peaks_cwp` function).
#' and (2) aligned (`obiAdjust_Rtime` function) peaks
#' @family dataset functions
#' @family retention time correction functions
#' @export
#' @examples
#' \dontrun{
#' file_name <-  system.file("extdata", "peakdet.rds", package = "AlpsLCMS")
#' peakdet <- base::readRDS(file_name)
#' opt_result_path <-  system.file("extdata", package = "AlpsLCMS")
#' preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#'
#' peakdet_align <- obiAdjust_Rtime(peakdet, params = preproc_params)
#' print(peakdet_align)
#' }


obiAdjust_Rtime <- function (peakdet, params) {

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

#' Peak Correspondence
#'
#' Peak correspondence is carried out by the *'lcms_groupChromPeaks'* method,
#' with parameters obtained form `IPO`. Peak Correspondece consist in
#' grouping peaks on retention time axis with the purspose of associate
#' them to spectra on the mass/chage axis.
#'
#' Note: signal processing generally  consists in three main steps:
#' (1) peak detection (`find_peaks_cwp` function),
#' (2) peak alignment (`obiAdjust_Rtime` function) and
#' (3) peak correspondence (`IPO_group_peaks` function).
#' The optimized set of parameters for this signal preprocessing can be obatained from `IPO` Package.
#'
#' After this stage the peak table is finally obtained.
#'
#' @param peakdet_align A lcms_dataset with (1) detected (`findChromPeaks_cwp`
#' function) and (2) aligned (`obiAdjust_Rtime` function) peaks.
#' @param params A converted parameters template from IPO parameters.
#' @return A lcms_dataset with (1) detected (`find_peaks_cwp` function),
#' (2) aligned (`align_rtime` function) and (3) grouped  (`IPO_group_peaks`
#' function) peaks.
#' @family dataset functions
#' @family peak correspondence functions
#' @export
#' @examples
#' \dontrun{
#' file_name <-  system.file("extdata", "peakdet_align.rds", package = "AlpsLCMS")
#' peakdet_align <- base::readRDS(file_name)
#' opt_result_path <-  system.file("extdata", package = "AlpsLCMS")
#' preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#'
#'
#' peak_table <- IPO_group_peaks(peakdet_align, params = preproc_params)
#' print(peak_table)
#' }
IPO_group_peaks <- function (peakdet_align, params) {
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

#' Filling missing values in a peak table
#'
#' In the imputation stage, the areas of the missing peaks of the peak table
#' that were not detected in the previous steps of the signal preprocessing workflow were integrated.
#'
#' @param peak_table A table of peaks with (possibly) missing values.
#' @return A peak table where the missing peaks have been filled.
#' @family dataset functions
#' @family imputation functions
#' @export
#' @examples
#' \dontrun{
#' file_name <-  system.file("extdata", "peak_table.rds", package = "AlpsLCMS")
#' peak_table <- base::readRDS(file_name)
#' peak_table_imputed <- lcms_fill_chrom_peaks(peak_table)
#'
#' print(peak_table_imputed)
#' }
lcms_fill_chrom_peaks <- function(peak_table){
  quiet <- function(x) {
    base::sink(base::tempfile())
    base::on.exit(base::sink())
    base::invisible(base::force(x))
  }

  cat("\n","Imputing missing peaks of the peak table.","\n")

  peak_table_imputed <-  base::suppressWarnings(
    base::suppressMessages(quiet(xcms::fillChromPeaks(peak_table))
    )
  )
  return(peak_table_imputed)
}



#' Total Ion Count (TIC)
#'
#' The function performs the Total Ion Count (TIC) for the polarity samples.
#' Function `lcms_tics` stores summarizes TIC information.
#'
#' NOTE: `lcms_tics` assumes that data is already filtered by polarity.
#'
#' @param dataset A lcms_dataset.
#' @param treatment Class groups of the samples.
#' @return Total Ion Count (TIC) for the polarity samples.
#' @family dataset functions
#' @family dataset_peak_table functions
#' @export
#' @examples
#' \dontrun{
#' dataset <- lcms_dataset_load(system.file
#'                                   ("extdata","dataset_pos.rds",
#'                                     package = "AlpsLCMS"))
#' tics <- lcms_tics(dataset)
#'
#' print(tics)
#' }
lcms_tics <- function(dataset, treatment = treatment){
  tics <- tibble::tibble(
    file = MSnbase::fromFile(dataset),
    fileName = Biobase::pData(dataset)$sampleNames[file],
    treatment = Biobase::pData(dataset)$treatment[file],
    ret_time = MSnbase::rtime(dataset),
    polarity = rep(unique(MSnbase::polarity(dataset),length(file))),
    tic = MSnbase::tic(dataset)
  )
  #Files sorted by treatment
  tics$fileName <- factor(tics$fileName,
                          levels = unique(tics$fileName)
                          [base::order(dataset$treatment)])
  return(tics)
}


#' Total Ion Count (TIC) plot
#'
#' The function performs a plot with the Total Ion Count (TIC).
#'
#' @param tics A Total Ion Count object generated with `lcms_tics`.
#' @param treatment Class groups of the samples.
#' @param rt Retention time boundaries.
#' @param plot_type The plot class, either boxplot or spectra.
#' @return Total Ion Count (TIC) for the polarity samples.
#' @family dataset functions
#' @family dataset_peak_table functions
#' @family chromatogram functions
#' @family visualization functions
#' @export
#' @examples
#' \dontrun{
#' dataset <- lcms_dataset_load(system.file
#'                                   ("extdata","dataset_pos.rds",
#'                                     package = "AlpsLCMS"))
#' tics <- lcms_tics(dataset)
#'
#' lcms_plot_tics(tics, treatment = treatment,
#'                rt = c(4, 8),plot_type = "spec")
#'
#' lcms_plot_tics(tics, treatment = treatment,
#'                rt = c(4, 8), plot_type = "boxplot")
#' }
lcms_plot_tics <- function(tics, treatment = treatment, rt = NULL, plot_type = "spec"){
  ret_time <- NULL
  tic <- NULL
  fileName <- NULL
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
    tics$fileName <- factor(sub(pattern = "(.*)\\..*$", replacement = "\\1", tics$fileName), levels = sub(pattern = "(.*)\\..*$", replacement = "\\1", levels(tics$fileName)))
    ggplot2::ggplot(tics) +
      ggplot2::geom_boxplot(ggplot2::aes(x = fileName, y = tic, fill = treatment)) +
      ggplot2::scale_fill_manual("Treatment", values = treatment_col) +
      ggplot2::scale_y_log10() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::ggtitle(paste("Boxplot of the Total Ion Count by sample", polarity))
  }
}


#' Retention Time Correction Plot
#'
#' It plots the retention time correction vs the original retention time for each of the samples
#' coloured by sample class.
#'
#' @param data An alignend lcms_dataset.
#' @return The plot for the retention time correction.
#' @family dataset functions
#' @family retention time correction functions
#' @family visualization functions
#' @export
#' @examples
#' \dontrun{
#' file_name <-  system.file("extdata", "peakdet_align.rds", package = "AlpsLCMS")
#' data <- base::readRDS(file_name)
#' rta_plot <- lcms_retention_time_alignment_plot(data)
#' print(rta_plot)
#' }
lcms_retention_time_alignment_plot <- function (data){
  fileIdx <- NULL
  ret_time_orig <- NULL
  ret_time_adj <- NULL
  treatment <- NULL
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

#' Base peak chromatogram
#'
#' Base peak chromatograms with retention time axis in minutes.
#'
#' @param chromatogram_object  A XChromatograms object.
#' @param treatment_col Color code by groups.
#' @param rtlim Retention time boundaries (e.g. c(4,8)).
#' @return A base peak chromatogram.
#' @family chromatogram functions
#' @family visualization functions
#' @export
#' @examples
#' \dontrun{
#' file_name <- system.file("extdata",
#'                          "dataset_pos.rds",
#'                           package = "AlpsLCMS")
#' dataset <- lcms_dataset_load(file_name)
#' treatment_col <- scales::hue_pal()(length(unique(dataset$treatment)))
#' names(treatment_col) <- unique(dataset$treatment)
#' base_peaks <- xcms::chromatogram(dataset, aggregationFun = "max")
#' lcms_plot_chrom(base_peaks, treatment_col, rtlim = c(4, 8))
#' }
lcms_plot_chrom <- function (chromatogram_object, treatment_col, rtlim = NULL){
  rtime <- NULL
  intensity  <- NULL
  min2sec <- 60
  message("Make sure that the column that contains the group class is called `treatment`")

  #we need to modify this in order to be more flexible (treatment_col)
  ret_times <- lapply(chromatogram_object, FUN = MSnbase::rtime)
  intensities <- lapply(chromatogram_object, FUN = MSnbase::intensity)

  graphics::plot(ret_times[[1]] / min2sec, intensities[[1]], type = "l",
       col = treatment_col[chromatogram_object$treatment][1], lwd = 1,
       xlab = "Retention time (min)", ylab = "Intensity  (A.U.)",
       xlim = rtlim,
       main = "Base Peak Chromatogram")

  for (i in 2:length(ret_times)){ # we need to modify the sequence using seq_along
    graphics::points(ret_times[[i]] / min2sec, intensities[[i]], lwd = 1,
           xlim = rtlim,
           type= "l", col = treatment_col[chromatogram_object$treatment][i])
    graphics::legend("topright", legend = names(treatment_col), fill = treatment_col)
  }
}

#' Image of Chromatographic Peaks by sample
#'
#' It plots the an image of the chromatographic peaks for each sample. This
#' function is useful if you are interested in knowing the effect of the
#' retention time correction on the chromatographic axis.
#'
#' @param x A lcms_dataset.
#' @param binSize Bin size for the the reteincion time axis.
#' @param xlim Retention time axis range. If Null, it uses the complete
#'   retention time range.
#' @param log If TRUE, it performs a logarithmic transformation on the data.
#' @param xlab Label for the retention time axis.
#' @param yaxt Options for the y axis. Use the default value.
#' @param main Plot title.
#' @param ... Additional parameters for image function.
#' @return An image with the detected chromatographic peak, for each sample.
#' @export
#' @family dataset functions
#' @family chromatogram functions
#' @family visualization functions
#' @examples
#' \dontrun{
#' file_name <-  system.file("extdata", "peakdet_align.rds", package = "AlpsLCMS")
#' dataset <- base::readRDS(file_name)
#' chr_peak_image <- lcms_plot_chrom_peak_image(dataset, binSize = 5, xlim = NULL, log = FALSE,
#'                                                xlab = "retention time (min)", yaxt = par("yaxt"))
#'                                                title(main ="Detected Peaks (Aligned)")
#' }

lcms_plot_chrom_peak_image<- function (x, binSize = 30, xlim = NULL, log = FALSE, xlab = "retention time",
                                       yaxt = graphics::par("yaxt"), main = "Chromatographic peak counts",
                                       ...)
{
  min2sec <- 60
  if (!methods::is(x, "XCMSnExp"))
    stop("'x' is supposed to be an 'XCMSnExp' object, but I got a ",
         class(x))
  if (is.null(xlim))
    xlim <- c(floor(min(MSnbase::rtime(x))), ceiling(max(MSnbase::rtime(x))))
  brks <- seq(xlim[1], xlim[2], by = binSize)
  if (brks[length(brks)] < xlim[2])
    brks <- c(brks, brks[length(brks)] + binSize)
  pks <- xcms::chromPeaks(x, rt = xlim)
  if (nrow(pks)) {
    rts <- split(pks[, "rt"], pks[, "sample"])
    cnts <- lapply(rts, function(z) {
      hst <- graphics::hist(z, breaks = brks, plot = FALSE)
      hst$counts
    })
    n_samples <- length(MSnbase::fileNames(x))
    sample_idxs <- 1:n_samples
    sample_idxs <- sample_idxs[!(as.character(sample_idxs) %in%
                                   names(rts))]
    if (length(sample_idxs)) {
      all_cnts <- vector("list", n_samples)
      all_cnts[as.numeric(names(cnts))] <- cnts
      zeros <- rep(0, (length(brks) - 1))
      all_cnts[sample_idxs] <- list(zeros)
      cnts <- all_cnts
    }
    cnts <- t(do.call(rbind, cnts))
    if (log)
      cnts <- log2(cnts)
    graphics::image(z = cnts, x = (brks - (brks[2] - brks[1])/2) / min2sec, xaxs = "r",
          xlab = xlab, yaxt = "n", ...)
    sample_labels <- stringr::str_extract(basename(MSnbase::fileNames(x)), "\\w+")
    graphics::axis(side = 2, at = seq(0, 1, length.out = n_samples),
         labels = FALSE)
    graphics::text(y = seq(0, 1, length.out = n_samples), graphics::par("usr")[1],
         cex = 0.6, labels = sample_labels,
         srt = 60, pos = 2, xpd = TRUE)
  }
  title(main = main)
}

#' Rearrange datafiles by class
#'
#' Although Data Preprocessing can be performed using any
#' filepath with the `xcms` Package, it is convenient to
#' rearrange data files by sample class (that is, all the
#' samples belonging to the same sample class are included
#' in the same folder). This is done because, after the
#' preprocessing stage, `MAIT` Package takes care
#' of Data Annotation and Metabolite Identification. `MAIT`
#' needs a specific directory structure for data managing.
#' Otherwise, it can't work properly.
#'
#' @param dataset A lcms_dataset.
#' @param dataDir A directory in where LC-MS files are
#' going to be saved.
#' @return LC-MS datafiles sorted by class treatment.
#' @family dataset functions
#' @family import/export functions
#' @export
#' @examples
#' \dontrun{
#' file_name <- system.file("extdata", "dataset_pos_rt_rs.rds", package = "AlpsLCMS")
#' dataset <- lcms_dataset_load(file_name)
#' path <- system.file("extdata","rearrange_mait", package = "AlpsLCMS")
#' lcms_rearrange_datafiles_by_class(dataset = dataset,
#'                             dataDir = path)
#' fileList = list.files(path = paste(path, list.files(path = path),
#'                            sep = "/"), full.names = TRUE)
#' print(fileList)
#' }
lcms_rearrange_datafiles_by_class <- function(dataset, dataDir) {
  Treatment <- NULL
  files <- Biobase::pData(dataset)$sampleNames
  files_treatment <- Biobase::pData(dataset)$treatment



  filetreat_info <- data.frame(Filename = files,
                               Treatment= files_treatment,
                               stringsAsFactors = FALSE)
  filetreat_info <- filetreat_info %>% dplyr::group_by(Treatment) %>% dplyr::arrange(Treatment)
  filetreat_info <- split(filetreat_info, filetreat_info$Treatment)

  drop_treatment <-function(x) {
    x$Treatment <- NULL
    colnames(x) <- "FileName"
    x
  }

  filetreat_info <- lapply(filetreat_info, FUN = drop_treatment)

  for (i in seq_along(filetreat_info)){
    filer <- filetreat_info[[i]][["FileName"]]
    foldname <-names(filetreat_info)[i]
    treatDir <- paste0(dataDir, "/", foldname, "/")
    if (length(dir(treatDir)) > 0) {
      cat("There are already directories / files in the folder. Not saving new ones. ")
      #cat("\n")
      return()
    }else{
      dir.create(treatDir)
      data_subset <- dataset %>% MSnbase::filterFile(file = filer)
      Biobase::fData(data_subset)$centroided <- TRUE
      Biobase::fData(data_subset)$peaksCount <- Biobase::fData(data_subset)$originalPeaksCount
      mzR::writeMSData(data_subset,
                       file = unlist(lapply(treatDir, FUN = paste0, filer)),
                       outformat = c("mzxml"), copy = FALSE)
    }
  }
}

#' Filter a mass window
#'
#' `filterMassAcquisition` is a function to filter a mass window when two
#' different mass windows have been acquired in the same experiment. Note that
#' the function works in a unique polarity, either positive or negative. Then,
#' the function extracts one of the sequential spectra within a `lcms_dataset`.
#' Spectra need to be alternating the "small and "large" mass window in the
#' file.
#'
#' @param dataset A lcms_dataset that contains two different mass ranges in two different acquisitions.
#' @param window Character. It refers to which mass range will be kept with
#'   `"small"` or `"large"`. `window = "small"` will keep the first spectrum of each
#'   couple of the intercalated spectra, while `window = "large"` will keep the
#'   second.
#' @return A `lcms_dataset` that only contains one mass range. The function
#'   filtered out the other spectra. For example, if "small" was set, The
#'   function will keep the first spectrum out of two because mass windows are
#'   intercalated within the file.
#' @export
#'
filterMassAcquisition <- function (dataset, window = "small"){
  if (window == "small"){
  number_Spec <- length(dataset@featureData@data[["polarity"]])
  selection = rep(1:2, times = number_Spec/2)
  Select <- dataset@featureData@data[["acquisitionNum"]]
  names(Select) <- selection
  small <- Select[names(Select)==1]
  large <- Select[names(Select)==2]
  dataset_window <- xcms::filterAcquisitionNum(dataset, small)
  }
  if (window == "large"){
    number_Spec <- length(dataset@featureData@data[["polarity"]])
    selection = rep(1:2, times = number_Spec/2)
    Select <- dataset@featureData@data[["acquisitionNum"]]
    names(Select) <- selection
    small <- Select[names(Select)==1]
    large <- Select[names(Select)==2]
    dataset_window <- xcms::filterAcquisitionNum(dataset, large)
  }
  return(dataset_window)
}

#' cromplot
#'
#' Function that creates four zoomed plots of a know metabolite to observe the
#' set the parameters for xcms preprocessing. Parameters such as peak width,
#' noise level, or potential neighbors can be studied.
#'
#' @param dataset a lcms_dataset
#' @param rt known retention time (minutes) of a metabolite that surely is
#'   within the dataset. Note that `rt` critically changes between columns.
#' @param mz known ionized mass of a metabolite that surely is within the dataset.
#' @param metabolite character. The name of the metabolite (for plotting).
#'
#' @return Plot of a known chromatographic peak
#' @export
#' @usage
#' cromplot(dataset,
#' rt = 5.8,
#' mz = 132.102,
#' metabolite = "(iso)leucine")
#'
cromplot <- function(dataset, rt, mz, metabolite = "metabolite"){
rt_max <- (rt*60)+120
rt_min <- (rt*60)-120
rtr <- c(rt_min, rt_max)

mz_max <- mz + 0.01
mz_min <- mz - 0.01
mzr <- c(mz_min, mz_max)

## extract the chromatogram
chr <- xcms::chromatogram(dataset, mz = mzr, rt = rtr)
MSnbase::plot(chr, col = "blue", main = metabolite)
MSnbase::plot(chr, col = "green", xlim = c((rt*60)-25, (rt*60)+25), ylim = c(0,1e+7), main = metabolite)
MSnbase::plot(chr, col = "orange", xlim = c((rt*60)-25, (rt*60)+25), ylim = c(5e+3, 1e+6), main = metabolite)
MSnbase::plot(chr, col = "red", xlim = c((rt*60)-25, (rt*60)+25), ylim = c(1e+3, 5e+5), main = metabolite)
}

#' Feature table
#'
#' `feature_values` is a wrapper from `xcms`. The function creates an intensity
#' matrix with all the features.
#'
#' @inheritParams xcms::featureValues
#' @inherit xcms::featureValues
#' @inheritSection xcms::featureValues @return
#' @export
#' @usage
#' feature_values(
#' object,
#' method = c("medret", "maxint", "sum"),
#' value = "into",
#' intensity = "into",
#' filled = TRUE,
#' missing = NA,
#' msLevel = integer()
#' )
feature_values <- function(...){
xdata <- xcms::featureValues(...)
xdata
}

#' Peak density parameters
#'
#' The method creates a template with the parameters for **correspondence**.
#' Extracted from the function `xcms::PeakDensityParam`.
#'
#' @inheritDotParams xcms::PeakDensityParam
#' @inherit xcms::PeakDensityParam
#'
#' @return PeakDensityParam S4 object
#' @export
#' @examples
#' p <- PeakDensityPar(sampleGroups = c(1,2))

PeakDensityPar <- function(...){
  object <- xcms::PeakDensityParam(...)
  object
}

#' Peak correspondence
#'
#' `groupPeaks` is a wrapper of the `xcms::groupChromPeaks` from `xcms`
#' package. It performs the **correspondence** (grouping) of chromatographic
#' peaks of the samples.
#'
#' @inheritDotParams xcms::groupChromPeaks
#' @inherit xcms::groupChromPeaks
#'
#' @return XCMSnExp object
#' @export
#'
groupPeaks <- function(...){
  object <- xcms::groupChromPeaks(...)
  object
}

#' Peak density parameters
#'
#' The method creates a template with the parameters for **alignment**.
#' Extracted from the function `xcms::PeakGroupsParam`.
#'
#' @inheritDotParams xcms::PeakGroupsParam
#'
#' @return PeakGroupsParam S4 object
#' @export
#'
PeakGroupsPar <- function(...){
  object <- xcms::PeakGroupsParam(...)
  object
}

#' Peak alignment
#'
#' `adjustRT` is a wrapper of the `xcms::adjustRtime` from `xcms` package. It
#' performs the **alignment** of chromatographic peaks of samples. There are two
#' methods that can be used depending on the parameters object.
#'
#' A) peakGroups: Correction of RT based on alignment of features (grouped
#' peaks) in most of samples. See [xcms::adjustRtime-peakGroups] for info.
#'
#' B) obiwarp: Alignment based on the complete mz-rt data. This method does not
#' require any identified peaks or defined features. See
#' [xcms::adjustRtime-obiwarp] for info.
#'
#' The parameter object [xcms::ObiwarpParam] or [xcms::PeakGroupsParam] created
#' for including in the `adjustRT` function defines which method will be used.
#'
#' @inheritDotParams xcms::adjustRtime
#' @inherit xcms::adjustRtime
#'
#' @return `XCMSnExp` object with peaks aligned. Note that `adjustRT` discard
#'   all previous peak grouping in the obtained object and thus, a new grouping
#'   process might be required.
#'   The [xcms::ObiwarpParam] function returns a
#'   `ObiwarpParam` class instance with all of the settings specified for
#'   obiwarp retention time adjustment and alignment.
#' @export
#'
adjustRT <- function(...){
  object <- xcms::adjustRtime(...)
  object
}

