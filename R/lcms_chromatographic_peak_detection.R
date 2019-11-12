#' Chromatographic peak detection (CentWave)
#'
#' The findChromPeaks_cwp function performs the chromatographic peak
#' detection on LC/GC-MS data. The standard method for peak detection
#' is *'CentWave'*. We must initialize its parameters according to the
#' `IPO` Package optimization. Peak detection aims to detect important
#' features (peaks) on the chromatographic axis. This will be useful
#' for a posterior peak alignment on the chormatophic axis.
#' Note: signal processing generally  consists in three main steps:
#' (1) peak detection (`lcms_findChromPeaks_cwp` function),
#' (2) peak alignment (`lcms_align_Rtime` function) and
#' (3) peak correspondence (`lcms_group_peaks` function).
#'
#' @param lcms_dataset An [lcms_dataset_family] object
#' @param params A converted parameters template from IPO parameters.
#' @examples
#' file_path <- system.file("extdata", "rearrange_mait", package = "NIHSlcms")
#' samples_mzxml <- list.files(file_path, recursive = TRUE, full.names = TRUE)
#' meta_path <- system.file("extdata", "metadata.xlsx", package = "NIHSlcms")
#' opt_result_path <-  system.file("extdata", package = "NIHSlcms")
#' lcms_preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#'
#' lcms_dataset <- suppressWarnings(lcms_read_samples(samples_mzxml, mode = "onDisk"))
#' metadata <- lcms_meta_read(meta_path)
#' lcms_dataset_meta <- lcms_meta_add(lcms_dataset, metadata, by = "sampleNames")
#'
#' peakdet <- lcms_findChromPeaks_cwp(lcms_dataset_meta, params = lcms_preproc_params)
#' print(peakdet)
#'
#' @return A lcms_dataset with detected peaks
#' @export
#'
lcms_findChromPeaks_cwp <- function (lcms_dataset, params) {
     quiet <- function(x) {
                base::sink(base::tempfile())
                base::on.exit(base::sink())
                base::invisible(base::force(x))
  }


  cat("\n","Finding chromatographic peaks using the optimized set of parameters obtained from IPO package.","\n")

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
                         quiet(xcms::findChromPeaks(lcms_dataset, param = cwp))
                                   )
                                  )
  return(peakdet)
}


