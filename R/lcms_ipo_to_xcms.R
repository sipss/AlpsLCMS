#' Read IPO parameters to XCMS format
#'
#' The function reads and converts from `IPO` to `XCMS` variable formats: The variable names in IPO and XCMS
#' presented some mismatches, thus we had to create a function called *lcms_read_ipo_to_xcms*
#' to achieve compatibility between packages when reading these parameters from a .csv file.
#'
#' @param opt_result_path A directory where the parameters file is stored
#'
#' @return A display of the chosen parameters
#' @export
#'
#' @examples
#' opt_result_path <-  system.file("extdata","ipo_opt", package = "NIHSlcms")
#' preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#' print(preproc_params)
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
#' presented some mismatches, thus we had to create a function called *lcms_convert_ipo_to_xcms*
#' to achieve compatibility between packages when these parameters are stored in a list.
#'
#' @param params A list with the IPO parameters, generated with the `lcms_write_opt_params` function
#' @return Parameters in XCMS format
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
