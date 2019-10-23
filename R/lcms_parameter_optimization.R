#' Default parameters for peak picking optimization
#'
#' The function creates default parameters for optimizing
#' peakpicking algorithms.We perform parameter optimization
#' on the XCMS preprocessing algorithms using the IPO Package.
#' This includes Peak Detection (‘Centwave’ and ‘Matched Filter’),
#' Retention Time Correction (‘obiwarp’) and Peak Correspondence
#' (‘Density’).
#' Use this function to generate the template within the
#' `lcms_peakpicking_optimization` function.
#'
#' @param noise numeric, minimum intensity needed to be included in the analysis.
#' @param snthresh numeric, set the signal to noise ratio threshold to be included.
#' @param min_peakwidth numeric(two values) with the expected minimal peak width in
#' chromatographic dimension. Set as a range (min, max) in seconds.
#' @param max_peakwidth numeric(two values) with the expected maximum peak width in
#' chromatographic dimension. Set as a range (min, max) in seconds.
#' @param optimize by default is TRUE. If FALSE, the function does not initialze the parameters
#' the parameters
#' @return A parameters template for peak picking optimization
#' @export
#' @family optimization functions
#'
#' @examples
#' default_peakpicking_params <- lcms_default_peakpicking_params(optimize = TRUE)
#' print(default_peakpicking_params)
#'
lcms_default_peakpicking_params <- function(noise = 5000, snthresh = 10,
                                            min_peakwidth = c(10, 30),
                                            max_peakwidth = c(35, 120),
                                            optimize = TRUE){
  if (optimize == TRUE){
    peakpickingParameters <- IPO::getDefaultXcmsSetStartingParams(method = c("centWave"))
    peakpickingParameters$noise <- noise
    peakpickingParameters$snthresh <- snthresh
    peakpickingParameters$min_peakwidth <- min_peakwidth
    peakpickingParameters$max_peakwidth <- max_peakwidth
  } else{
    peakpickingParameters <- NULL
  }
  peakpickingParameters
}

#' Peak picking optimization
#'
#' The function optimize parameters considering a set of samples
#' for the peak picking algorithm using the IPO Package.
#' This includes Peak Detection (‘Centwave’ and ‘Matched Filter’),
#' Retention Time Correction (‘obiwarp’) and Peak Correspondence
#' (‘Density’).
#'
#' @param lcms_dataset A [lcms_dataset_family] object
#' @param peakpickingParameters Parameters for peak picking
#' @param opt_path Path where optimization samples and optimized parameters are save. If NULL, a temporary folder is created.
#' @param nSlaves Number of slaves the optimization process should spawn.
#' @param plots Defines if plots should be generated (TRUE) or not (FALSE) in a subfolder called "plot_ipo".
#' @return A peak picking list with the best setting
#' @export
#' @family optimization functions
#' @examples
#' \dontrun{
#' file_name <- system.file("extdata", "lcms_dataset_rt_pos_rs.rds", package = "NIHSlcms")
#' lcms_dataset <- lcms_dataset_load(file_name)
#' default_peakpicking_params <- lcms_default_peakpicking_params(optimize = TRUE)
#' resultPeakpicking <- lcms_peakpicking_optimization(lcms_dataset,
#'                                                    default_peakpicking_params,
#'                                                    opt_path = NULL)
#' print(resultPeakpicking)}
#'
lcms_peakpicking_optimization <- function (lcms_dataset, peakpickingParameters, nSlaves = 1, opt_path, plots = TRUE){
  message("This function requires a folder without previous LC-MS data")

  if(is.null(peakpickingParameters)){
    resultPeakpicking <- NULL
  } else{
    filenames <- Biobase::pData(lcms_dataset)$sampleNames
    filer <- filenames
    former_dir <- getwd()
    if(is.null(opt_path)){
      opt_path <-tempdir()

    }
    setwd(opt_path)

    ## Get the spectra
    data_subset <- lcms_dataset %>% MSnbase::filterFile(file = filenames)
    Biobase::fData(data_subset)$centroided <- TRUE
    Biobase::fData(data_subset)$peaksCount <- Biobase::fData(data_subset)$originalPeaksCount
    print("Be aware: do not run twice using the same output directory")
    print("The algorithm is not able to rewrite files that are already in the directory")

    mzR::writeMSData(data_subset, file = filer, outformat = c("mzxml"), copy = FALSE)

    print("Saving filtered chromatogram...")

    #samples_op <- fs::dir_ls(opt_path , glob = "*.mzXML")
    print("Performing peak detection parameter optimization. This will take some time...")
    time.xcmsSet <- system.time({ # measuring time
      base::suppressWarnings(
        base::suppressMessages(
          resultPeakpicking <- IPO::optimizeXcmsSet(files =  opt_path,
                                               params = peakpickingParameters,
                                               nSlaves = nSlaves,
                                               subdir = "plot_ipo",
                                               plot = plots)
        )
      )
    })
    setwd(former_dir)
  }
  return(resultPeakpicking)
}

#' Default parameters for optimization of retention time correction and grouping parameters
#'
#' The function creates default parameters for optimizing
#' retention time correction and grouping algorithms.
#' We perform parameter optimization on the XCMS preprocessing
#' algorithms using the IPO Package.
#' Use this function to create the template within `lcms_corgroup_optimization` function.
#'
#' @param profStep set the m/z step for generating profile (matrix) data from raw mass spectral data.
#' Smaller steps yield more precision at the cost of greater memory usage.
#' @param gapExtend numeric(1) defining the penalty for gap enlargement.
#' The default value for gapExtend depends on the value of distFun, for distFun = "cor"
#' and distFun = "cor_opt" it is 2.4, for distFun = "cov" 11.7, for distFun = "euc" 1.8
#' and for distFun = "prd" 7.8.
#' @param optimize by default is TRUE. If FALSE, the function does not optimize
#' the parameters
#' @return A parameters template for retention time correction and grouping optimization
#' @export
#' @family optimization functions
lcms_default_corgroup_params <- function(profStep = 1, gapExtend = 2.7, optimize = TRUE){

  if (optimize == TRUE){
    retcorGroupParameters <- IPO::getDefaultRetGroupStartingParams()
    retcorGroupParameters$profStep <- profStep
    retcorGroupParameters$gapExtend <- gapExtend
  } else{
    retcorGroupParameters <- NULL
  }
  retcorGroupParameters
}

#' Optimization of retention time correction and grouping parameters.
#'
#' The function optimize parameters considering a set of samples
#' for the retention time correction and grouping using the IPO Package.
#'
#' @param optimizedXcmsSetObject XCMS object conatining the `best_settings` parameters.
#' This object may be created after running `lcms_peakpicking_optimization` and extract
#' best_settings$xset (e.g. optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset)
#' @param retcorGroupParameters Parameters for retention time correction and optimization
#' @param nSlaves Number of slaves the optimization process should spawn.
#' @param plots Defines if plots should be generated (TRUE) or not (FALSE) in a
#' subfolder called "plot_ipo".
#' @return a list with the optimization of parameters for retention time and grouping.
#' @export
#' @family optimization functions.
lcms_corgroup_optimization <- function (optimizedXcmsSetObject,
                                        retcorGroupParameters,
                                        nSlaves = 1,
                                        opt_path,
                                        plots = TRUE){

  if(is.null(optimizedXcmsSetObject) | is.null(retcorGroupParameters)){
    resultRetcorGroup <- NULL
  } else{
    former_dir <- getwd()
    setwd(opt_path)
    print("Performing retention time and grouping
        parameter optimization. This will take some time...")
    time.RetGroup <- system.time({ # measuring time
      base::suppressWarnings(
        base::suppressMessages(
          resultRetcorGroup <-
            IPO::optimizeRetGroup(xset = optimizedXcmsSetObject,
                             params = retcorGroupParameters,
                             nSlaves = nSlaves,
                             subdir = "plot_ipo",
                             plot = plots)
        )
      )
    })
    setwd(former_dir)
  }
  return(resultRetcorGroup)
}

#' Displaying and Storing optimized settings
#'
#' The function allows visulizing the parameter optimization results by `IPO` Package
#' allows in the RStudio console. Also you can save this results in plain text files
#' (i.e. a .CVS file).
#' @param results_pp object from the `lcms_peakpicking_optimization`function
#' @param results_rtcg object from the `lcms_corgroup_optimization`function
#' @param opt_result_path A directory to save the parameters file
#' @param csv if TRUE, it writes a file in csv format
#' @param console if TRUE, it displays the params on the console
#'
#' @return A file with the params from the IPO optimization
#' @export
#' @family optimization functions
#' @examples
#' \dontrun{
#' opt_result_path <- "C:/my_directory"
#' write_opt_params(resultPeakpicking, resultRetcorGroup, opt_result_path)
#' }
write_opt_params<- function(results_pp,
                            results_rtcg,
                            opt_result_path,
                            csv = TRUE,
                            console = TRUE){
  if (is.null(results_pp) | is.null(results_rtcg)){
    paramsPP <- list()
    paramsPP$min_peakwidth <- 20
    paramsPP$max_peakwidth <- 50
    paramsPP$ppm <- 25
    paramsPP$mzdiff <- -0.001
    paramsPP$snthresh <- 10
    paramsPP$noise <- 0
    paramsPP$prefilter <- 3
    paramsPP$value_of_prefilter <- 100
    paramsPP$mzCenterFun <-"wMean"
    paramsPP$integrate <- 1
    paramsPP$fitgauss <-FALSE
    paramsPP$verbose.columns <- FALSE

    paramsRTCGroup <- list()
    paramsRTCGroup$retcorMethod = "obiwarp"
    paramsRTCGroup$plottype = "none"
    paramsRTCGroup$profStep <- 1
    paramsRTCGroup$center <-  NULL
    paramsRTCGroup$response <- 1
    paramsRTCGroup$distFunc <- "cor_opt"
    paramsRTCGroup$gapInit <- NULL
    paramsRTCGroup$gapExtend <- NULL
    paramsRTCGroup$factorDiag <- 2
    paramsRTCGroup$factorGap <- 1
    paramsRTCGroup$localAlignment <- 0
    paramsRTCGroup$initPenalty <- 0

    paramsRTCGroup$bw <- 30
    paramsRTCGroup$minfrac <- 0.5
    paramsRTCGroup$minsamp <- 1
    paramsRTCGroup$max <- 50
    paramsRTCGroup$mzwid <- 0.25

  } else{
    paramsPP <- results_pp$best_settings$parameters
    paramsRTCGroup <- results_rtcg$best_settings

  }


  if (console == TRUE){
    IPO::writeRScript(paramsPP, paramsRTCGroup)
  }
  if (csv == TRUE){
    IPO::writeParamsTable(paramsPP, paramsRTCGroup,
                     sep = ",", paste0(opt_result_path, "/params.csv"))
  }
  paramsPP
}
