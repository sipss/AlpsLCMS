#' Default parameters for peak picking optimization
#'
#' The function creates default parameters for optimizing
#' peakpicking algorithms.Parameter optimization
#' on the XCMS preprocessing algorithms was performed using the IPO Package.
#' This includes Peak Detection (‘Centwave’),
#' Retention Time Correction (‘obiwarp’) and Peak Correspondence
#' (‘Density’). Use this function to generate the template within the
#' `lcms_peakpicking_optimization` function.
#'
#' @param noise numeric, minimum intensity needed to be included in the analysis.
#' @param snthresh numeric, set the signal to noise ratio threshold to be included.
#' @param min_peakwidth numeric (two values) with the expected minimal peak width in
#' chromatographic dimension. Set as a range (min, max) in seconds.
#' @param max_peakwidth numeric (two values) with the expected maximum peak width in
#' chromatographic dimension. Set as a range (min, max) in seconds.
#' @param optimize by default is TRUE. If FALSE, the function does not initialze the parameters
#' the parameters
#' @return A parameters template for peak picking optimization
#' @family optimization functions
#' @export
#' @examples
#' \dontrun{
#' default_peakpicking_params <- lcms_default_peakpicking_params(optimize = TRUE)
#' print(default_peakpicking_params)
#' }
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
#' This includes Peak Detection (‘Centwave),
#' Retention Time Correction (‘obiwarp’) and Peak Correspondence
#' (‘Density’).
#'
#' @param dataset A lcms_dataset
#' @param peakpickingParameters Parameters for peak picking
#' @param opt_path Path where optimization samples are saved.
#' @param nSlaves Number of slaves the optimization process should spawn.
#' @param plots Defines if plots should be generated (TRUE) or not (FALSE) in a subfolder called "plot_ipo" (default).
#' @param subdir Folder where surface plots are saved. If NULL they are displayed by the graphical device.
#' @return A peak picking list with the best setting.
#' @family dataset functions
#' @family optimization functions
#' @family visualization functions
#' @export
#' @examples
#' \dontrun{
#' opt_path <- system.file("extdata", "ipo_opt", package = "NIHSlcms")
#' file_name <- system.file("extdata", "dataset_pos_rt_rs.rds", package = "NIHSlcms")
#' dataset <- lcms_dataset_load(file_name)
#' default_peakpicking_params <- lcms_default_peakpicking_params(optimize = TRUE)
#' result_peakpicking <- lcms_peakpicking_optimization(dataset,
#'                                                    default_peakpicking_params,
#'                                                    opt_path = opt_path,
#'                                                    subdir = NULL)
#' }
lcms_peakpicking_optimization <- function (dataset, peakpickingParameters,
                                           nSlaves = 1, opt_path, subdir ="plot_ipo",
                                           plots = TRUE){

  if(is.null(peakpickingParameters)){
    resultPeakpicking <- NULL
  } else{
    filenames <- Biobase::pData(dataset)$sampleNames
    filer <- filenames
    former_dir <- getwd()
    setwd(opt_path)

    ## Get the spectra
    data_subset <- dataset %>% MSnbase::filterFile(file = filenames)
    Biobase::fData(data_subset)$centroided <- TRUE
    Biobase::fData(data_subset)$peaksCount <- Biobase::fData(data_subset)$originalPeaksCount
    cat("Be careful if you run twice the function using the same output directory.", "\n")
    cat("The algorithm will remove the files stored in it, and then write the new files.", "\n")

    mzxml_in_opt_path <- lcms_list_mzxml_samples(opt_path, file_format = "mzXML",
                                                 rawconverter_path = NULL)

    file_names_opt_path <- stringr::str_c(stringr::str_match(mzxml_in_opt_path,"\\w+\\.mzXML$"))
    num_mzxml_opt_path <- length(file_names_opt_path)
    file_names_union <- union(file_names_opt_path, filer)

    if(num_mzxml_opt_path > 0){
      file.remove(file_names_opt_path)
    }

    mzR::writeMSData(data_subset, file = filer, outformat = c("mzxml"), copy = FALSE)
    aux_filer <- stringr::str_c(filer,collapse = " ")
    cat(stringr::str_c("Samples used for optimization:",
                       "\n", "\t",aux_filer, "\n",collapse =" "))

    cat("Saving filtered chromatogram...","\n")

    #samples_op <- fs::dir_ls(opt_path , glob = "*.mzXML")
    cat("Performing retention time and grouping
        parameter optimization. This will take some time...","\n")
    time.xcmsSet <- system.time({ # measuring time
      base::suppressWarnings(
        base::suppressMessages(
          resultPeakpicking <- IPO::optimizeXcmsSet(files =  opt_path,
                                               params = peakpickingParameters,
                                               nSlaves = nSlaves,
                                               subdir = subdir,
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
#' Parameter optimization on the XCMS preprocessing
#' algorithms is performed using the IPO Package.
#' Use this function to create the template within `lcms_corgroup_optimization` function.
#'
#' @param profStep set the m/z step for generating profile (matrix) data from raw mass spectral data.
#' Smaller steps yield more precision at the cost of greater memory usage.
#' @param gapExtend numeric(1) defining the penalty for gap enlargement.
#' The default value for gapExtend depends on the value of distFun, for distFun = "cor"
#' and distFun = "cor_opt" it is 2.4, for distFun = "cov" 11.7, for distFun = "euc" 1.8
#' and for distFun = "prd" 7.8.
#' @param optimize by default is TRUE. If FALSE, the function does not optimize
#' the parameters.
#' @return A parameters template for retention time correction and grouping optimization
#' @family optimization functions.
#' @export
#' @examples
#' \dontrun{
#' default_retcorgroup_params <- lcms_default_retcorgroup_params(optimize = TRUE)
#' print(default_retcorgroup_params)
#' }
lcms_default_retcorgroup_params <- function(profStep = 1,
                                            gapExtend = 2.7,
                                            optimize = TRUE){

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
#' The function optimizes parameters considering a set of samples
#' for the retention time correction and grouping using the IPO Package.
#'
#' @param optimizedXcmsSetObject XCMS object conatining the `best_settings` parameters.
#' This object may be created after running `lcms_peakpicking_optimization` and extract
#' best_settings$xset (e.g. optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset)
#' @param retcorGroupParameters Parameters for retention time correction and optimization
#' @param nSlaves Number of slaves the optimization process should spawn.
#' @param opt_path Path where optimization samples are saved.
#' @param subdir Folder where surface plots are save. If NULL they are displayed by the graphical device.
#' @param plots Defines if plots should be generated (TRUE) or not (FALSE) in a
#' subfolder called "plot_ipo"(default).
#' @return A list with the optimization of parameters for retention time and grouping.
#' @family optimization functions
#' @family visualization functions
#' @export
#' @examples
#' \dontrun{
#' file_name_pp <- system.file("extdata", "result_peakpicking.rds", package = "NIHSlcms")
#' optimizedXcmsSetObject <-base::readRDS(file_name_pp)$best_settings$xset
#' file_name_rcg <- system.file("extdata", "default_retcorgroup_params.rds", package = "NIHSlcms")
#' default_retcorgroup_params <- base::readRDS(file_name_rcg)
#' opt_path <-  system.file("extdata", package = "NIHSlcms")
#'
#' result_retcorgroup <- lcms_retcorgroup_optimization(optimizedXcmsSetObject,
#'                                                    default_retcorgroup_params,
#'                                                    opt_path = opt_path,
#'                                                    subdir = NULL)
#' }
lcms_retcorgroup_optimization <- function (optimizedXcmsSetObject,
                                        retcorGroupParameters,
                                        nSlaves = 1,
                                        opt_path,
                                        subdir ="plot_ipo",
                                        plots = TRUE){

  if(is.null(optimizedXcmsSetObject) | is.null(retcorGroupParameters)){
    resultRetcorGroup <- NULL
  } else{
    former_dir <- getwd()
    setwd(opt_path)

    quiet <- function(x) {
           base::sink(base::tempfile())
           base::on.exit(base::sink())
           base::invisible(base::force(x))
        }


    cat("Performing retention time and grouping
        parameter optimization. This will take some time...","\n")

   time.RetGroup <- system.time({ # measuring time
      base::suppressWarnings(
        base::suppressMessages(
          resultRetcorGroup <-quiet(
                    IPO::optimizeRetGroup(xset = optimizedXcmsSetObject,
                                          params = retcorGroupParameters,
                                          nSlaves = nSlaves,
                                          subdir = subdir,
                                          plot = plots)
        )
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
#'
#' @param results_pp object from the `lcms_peakpicking_optimization` function. If NULL, the default parameters for xcms are loaded.
#' @param results_rtcg object from the `lcms_corgroup_optimization` function.
#' @param opt_result_path A directory to save the parameters file.
#' @param csv if TRUE, it writes a file in csv format.
#' @param console if TRUE, it displays the params on the console.
#' @return A file with the params from the IPO optimization.
#' @family optimization functions
#' @family import/export functions
#' @export
#' @examples
#' \dontrun{
#'  opt_result_path <-system.file("extdata", "ipo_opt", package = "NIHSlcms")
#'  file_name_pp <- system.file("extdata", "result_peakpicking.rds", package = "NIHSlcms")
#'  file_name_rcg <- system.file("extdata", "result_retcorgroup.rds", package = "NIHSlcms")
#'  result_peakpicking <- base::readRDS(file_name_pp)
#'  result_retcorgroup <- base::readRDS(file_name_rcg)
#'
#'  lcms_write_opt_params(result_peakpicking, result_retcorgroup, opt_result_path)
#'  }
lcms_write_opt_params<- function(results_pp,
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
}

#' Read IPO parameters to XCMS format
#'
#' The function reads and converts from `IPO` to `XCMS` variable formats: The variable names in IPO and XCMS
#' presented some mismatches, the function *lcms_read_ipo_to_xcms* was created
#' to achieve compatibility between packages when reading these parameters from a .csv file.
#'
#' @param opt_result_path A directory where the parameters file is stored.
#' @return A display of the chosen parameters.
#' @family optimization functions
#' @family import/export functions
#' @export
#' @examples
#' \dontrun{
#' opt_result_path <-  system.file("extdata","ipo_opt", package = "NIHSlcms")
#' preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#' print(preproc_params)
#' }
lcms_read_ipo_to_xcms <- function(opt_result_path){
  params <- utils::read.csv(paste0(opt_result_path, "/params.csv"), stringsAsFactors = FALSE)
  preproc_params <- list(ppm = params$ppm,
                         peakwidth = c(params$min_peakwidth, params$max_peakwidth),
                         snthresh = as.numeric(params$snthresh),
                         prefilter = c(params$prefilter, params$value_of_prefilter),
                         mzCenterFun = params$mzCenterFun,
                         integrate =  params$integrate,
                         mzdiff = params$mzdiff,
                         fitgauss = params$fitgauss,
                         noise = params$noise,
                         verbose.columns = params$verbose.columns,
                         profStep = as.numeric(params$profStep),
                         centerSample = params$center,
                         response = as.numeric(params$response),
                         distFun = params$distFunc,
                         gapInit = as.numeric(params$gapInit),
                         gapExtend = as.numeric(params$gapExtend),
                         factorDiag = as.numeric(params$factorDiag),
                         factorGap = as.numeric(params$factorGap),
                         localAlignment = as.logical(params$localAlignment),
                         initPenalty = 0,
                         bw = as.numeric(params$bw),
                         minFraction = as.numeric(params$minfrac),
                         minSamples = as.numeric(params$minsamp),
                         maxFeatures = as.numeric(params$max),
                         mzwid = as.numeric(params$mzwid))
  preproc_params

}

#' Converts IPO parameters to XCMS format
#'
#' The function converts from `IPO` to `XCMS` variable formats: The variable names in IPO and XCMS
#' presented some mismatches, thus the function *lcms_convert_ipo_to_xcms* was created
#' to achieve compatibility between packages when these parameters are stored in a list.
#'
#' @param params A list with the IPO parameters, generated with the `lcms_write_opt_params` function.
#' @return Parameters in XCMS format
#' @family optimization functions
#' @family import/export functions
#' @export
#' @examples
#' \dontrun{
#' opt_result_path <-  system.file("extdata","ipo_opt", "params.rds", package = "NIHSlcms")
#' params <- base::readRDS(opt_result_path)
#' preproc_params <- lcms_convert_ipo_to_xcms(params)
#' print(preproc_params)
#' }

lcms_convert_ipo_to_xcms <- function(params){
  preproc_params <- list(ppm = params$ppm,
                         peakwidth = c(params$min_peakwidth, params$max_peakwidth),
                         snthresh = as.numeric(params$snthresh),
                         prefilter = c(params$prefilter, params$value_of_prefilter),
                         mzCenterFun = params$mzCenterFun,
                         integrate =  params$integrate,
                         mzdiff = params$mzdiff,
                         fitgauss = params$fitgauss,
                         noise = params$noise,
                         verbose.columns = params$verbose.columns,
                         profStep = as.numeric(params$profStep),
                         centerSample = params$center,
                         response = as.numeric(params$response),
                         distFun = params$distFunc,
                         gapInit = as.numeric(params$gapInit),
                         gapExtend = as.numeric(params$gapExtend),
                         factorDiag = as.numeric(params$factorDiag),
                         factorGap = as.numeric(params$factorGap),
                         localAlignment = as.logical(params$localAlignment),
                         initPenalty = 0,
                         bw = as.numeric(params$bw),
                         minFraction = as.numeric(params$minfrac),
                         minSamples = as.numeric(params$minsamp),
                         maxFeatures = as.numeric(params$max),
                         mzwid = as.numeric(params$mzwid))
  preproc_params

}
