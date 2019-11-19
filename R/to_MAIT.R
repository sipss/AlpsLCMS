#' Write Parameters Table for MAIT functionalities
#'
#' @param listParameters Parameters
#' @param folder A directory to store a .csv file where the  parameter table is stored.
#' @return A template with MAIT parameters.
#' @export
#' @examples
#' \dontrun{
#' dataDir <-  system.file("extdata", "rearrange_mait", package = "NIHSlcms")
#' MAITparams_path <- system.file("extdata", "Results_project", package = "NIHSlcms")
#' opt_result_path <-  system.file("extdata", package = "NIHSlcms")
#' lcms_preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#' parameters <- c(dataDir = dataDir, lcms_preproc_params)
#'
#' quiet <- function(x) {
#'              base::sink(base::tempfile())
#'              base::on.exit(base::sink())
#'              base::invisible(base::force(x))}
#'
#' quiet(getClassDef("MAIT","MAIT"))
#' MAIT.object = signature(MAIT.Object = "MAIT")
#' MAIT.object <- methods::new("MAIT")
#' MAIT.object@RawData@parameters@sampleProcessing <- parameters
#' lcms_writeParameterTable(MAIT::parameters(MAIT.object), folder = MAITparams_path)
#' }
#'
lcms_writeParameterTable <- function(listParameters, folder){
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
#' `MAIT` Package allows to annotate the peaks of the peak table provided by `XCMS`.
#' This is done in three different stages:
#'
#' * Annotation with `CAMERA`
#' * Annotation using predefined biotransformation
#' * Annotation using the Human Metabolome Database (HMDB)
#'
#' *Note: The dataset must be converted from an object of the XCMS Package
#' to an object of the MAIT Package.
#'
#' @param dataDir directory with LC-MS datafiles
#' @param projectDir path to the the project directory
#' @param project name of the
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
#' projectDir <-  system.file("extdata", package = "NIHSlcms")
#' project <- "project"
#'
#' peak_table_MAIT <- lcms_to_MAIT(dataDir = dataDir,
#'                                 projectDir = projectDir,
#'                                 project = project,
#'                                 preproc_params = lcms_preproc_params,
#'                                 peakTable = peak_table)
#' print(peak_table_MAIT)
lcms_to_MAIT <- function (dataDir = NULL, projectDir = NULL, project = NULL, preproc_params = NULL,  peakTable = NULL){

  quiet <- function(x) {
    base::sink(base::tempfile())
    base::on.exit(base::sink())
    base::invisible(base::force(x))
  }

  cat("Building a MAIT object from xcms data.","\n")
  parameters <- c(dataDir = dataDir, preproc_params)
  getClassDef("MAIT","MAIT")
  MAIT.object = signature(MAIT.Object = "MAIT")
  MAIT.object <- methods::new("MAIT")




  if (is.null(dataDir)) {
    stop("No input directory was given")
  }
  if (is.null(preproc_params)) {
    stop("No parameters for preprocessing were included")
  }
  if (is.null(peakTable)) {
    stop("No peak was included")
  }
  print(MAIT::resultsPath(MAIT.object))
  MAIT.object@RawData@parameters@sampleProcessing <- parameters


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

  if (!is.null(projectDir)) {
    resultsPath <- paste(projectDir, paste("Results", project, sep = "_"), sep ="/")
    if (any(base::dir.exists(resultsPath))) {
      cat("There are already directories / files in the folder. Not saving new ones.")
      cat("\n")
      }else{
        dir.create(resultsPath)
        peakTable <- base::suppressMessages(as(peakTable, "xcmsSet"))
        peakTable <-  base::suppressMessages(
          base::suppressWarnings(quiet(xcms::fillPeaks(peakTable))
          )
        )
      }
  }else{
        resultsPath <- "Results"
        dir.create(resultsPath)
  }

  lcms_writeParameterTable(base::suppressMessages(base::suppressWarnings(
                                                                        quiet(MAIT::parameters(MAIT.object)))),
                          folder = resultsPath)

  peakTable <- base::suppressMessages(as(peakTable, "xcmsSet"))
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
