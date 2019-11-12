#' Write Parameters Table for MAIT functionalities
#'
#'`MAIT` Package allows to annotate the peaks of the peak table provided by `XCMS`.
#'This is done in three different stages:
#' * Annotation with `CAMERA`
#' * Annotation using predefined biotransformation
#' * Annotation using the Human Metabolome Database (HMDB)
#'
#' *Note: The dataset must be converted from an object of the XCMS Package
#' to an object of the MAIT Package.
#'
#' @param listParameters Parameters
#' @param folder A directory
#' @return A template with MAIT parameters. Also write the `MAITparameters.csv´ file.
#' @export
#'
writeParameterTable <- function(listParameters, folder){
  outputTable <- as.matrix(c(unlist(listParameters@sampleProcessing),
                             unlist(listParameters@peakAnnotation),
                             unlist(listParameters@peakAggregation),
                             unlist(listParameters@sigFeatures),
                             unlist(listParameters@biotransformations),
                             unlist(listParameters@identifyMetabolites),
                             unlist(listParameters@classification),
                             unlist(listParameters@plotPCA),
                             unlist(listParameters@plotHeatmap)))
  colnames(outputTable) <- c("Value")
  write.csv(file=paste(folder,"MAITparameters.csv", sep="/"), outputTable)
  return(outputTable)
}

#' to MAIT
#'
#' Function to create a MAIT object using data from xcms
#'
#' @param dataDir directory with LC-MS datafiles
#' @param project name of the project
#' @param preproc_params params from XCMS
#' @param peakTable XCMS procesed LC-MS data
#'
#' @return A MAIT object
#' @export
#' @examples
#'
#' file_name <-  system.file("extdata", "peak_table_imputed.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' opt_result_path <-  system.file("extdata", package = "NIHSlcms")
#' lcms_preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#'
#' dataDir <-  system.file("extdata", "rearrange_mait", package = "NIHSlcms")
#' project <- "project"
#'
#' peak_table_MAIT <- lcms_to_MAIT(dataDir = dataDir,
#'                                 project = project,
#'                                 preproc_params = lcms_preproc_params,
#'                                 peakTable = peak_table)
#' print(peak_table_MAIT)
lcms_to_MAIT <- function (dataDir = NULL, project = NULL, preproc_params = NULL,  peakTable = NULL){

  quiet <- function(x) {
    base::sink(base::tempfile())
    base::on.exit(base::sink())
    base::invisible(base::force(x))
  }

  cat("\n","Building a MAIT object from xcms data.","\n")
  parameters <- c(dataDir = dataDir, preproc_params)
  getClassDef("MAIT","MAIT")
  MAIT.object = signature(MAIT.Object = "MAIT")
  MAIT.object <- methods::new("MAIT")

  if (is.null(dataDir)) {
    stop("No input directory was given")
  }
  if (is.null(project)) {
    stop("No project name was included")
  }
  if (is.null(preproc_params)) {
    stop("No parameters for preprocessing were included")
  }
  if (is.null(peakTable)) {
    stop("No peak was included")
  }


  MAIT.object@RawData@parameters@sampleProcessing <- parameters
  NIHSlcms::writeParameterTable(base::suppressMessages(base::suppressWarnings(
                                                              quiet(MAIT::parameters(MAIT.object)))),
                                                              folder = base::suppressMessages(
                                                                                base::suppressWarnings(quiet(MAIT::resultsPath(MAIT.object))))) # warning: resultsPath
  class <- list.files(dataDir)
  classNum <- vector(length = length(class))
  fileList <- list.files(path = paste(dataDir, list.files(path = dataDir),
                                      sep = "/"), full.names = TRUE)
  for (i in 1:length(class)) {
    classNum[i] <- length(list.files(paste(dataDir, list.files(dataDir)[i],
                                           sep = "/")))
  }

  classes <- rep(class, classNum)

  if (length(list.files(dataDir)) == 1) {
    warning("Warning: Input data only has one class!")
  }
  if (is.null(project)) {
    warning("Warning: Project name is empty!")
  }
  if (!is.null(project)) {
    resultsPath <- paste("Results", project, sep = "_")
    dir.create(resultsPath)
  }
  else {
    resultsPath <- "Results"
    dir.create(resultsPath)
  }

  peakTable <- as(peakTable, "xcmsSet")
  peakTable <-  base::suppressMessages(
                      base::suppressWarnings(quiet(xcms::fillPeaks(peakTable))
                                            )
                                   )

  fPeaks <- list(peakTable)
  names(fPeaks) <- "xcmsSet"
  MAIT.object@RawData@data <- fPeaks
  MAIT.object@PhenoData@classes <- class
  MAIT.object@PhenoData@classNum <- classNum
  MAIT.object@PhenoData@resultsPath <- resultsPath
  MAIT <- MAIT.object
}
