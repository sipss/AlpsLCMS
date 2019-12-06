#' Write Parameters Table for MAIT functionalities
#'
#' @param listParameters Parameters
#' @param folder A directory to store a .csv file where the  parameter table is stored.
#' @return A template with MAIT parameters.
#' @export
#' @examples
#' \dontrun{
#' data_dir <-  system.file("extdata", "rearrange_mait", package = "NIHSlcms")
#' mait_params_path <- system.file("extdata", "results_project", package = "NIHSlcms")
#' opt_result_path <-  system.file("extdata", package = "NIHSlcms")
#' preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#' parameters <- c(dataDir = dataDir, preproc_params)
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
#' lcms_write_parameter_table(MAIT::parameters(MAIT.object), folder = mait_params_path)
#' }
lcms_write_parameter_table <- function(listParameters, folder){
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
#' @param data_dir directory with LC-MS datafiles
#' @param project_dir path to the the project directory
#' @param project name of the
#' @param preproc_params params from XCMS
#' @param peak_table XCMS procesed LC-MS data
#'
#' @return A MAIT object
#' @export
#' @examples
#' file_name <-  system.file("extdata", "peak_table_imputed.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' opt_result_path <-  system.file("extdata", package = "NIHSlcms")
#' preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#'
#' data_dir <-  system.file("extdata", "rearrange_mait", package = "NIHSlcms")
#' project_dir <-  system.file("extdata", package = "NIHSlcms")
#' project <- "project"
#'
#' peak_table_mait <- lcms_to_mait(data_dir = data_dir,
#'                                 project_dir = project_dir,
#'                                 project = project,
#'                                 preproc_params = preproc_params,
#'                                 peak_table = peak_table)
#' print(peak_table_mait)
#'
lcms_to_mait <- function (data_dir = NULL, project_dir = NULL, project = NULL, preproc_params = NULL,  peak_table = NULL){

  quiet <- function(x) {
    base::sink(base::tempfile())
    base::on.exit(base::sink())
    base::invisible(base::force(x))
  }

  cat("Building a MAIT object from xcms data.","\n")
  parameters <- c(dataDir = data_dir, preproc_params)
  getClassDef("MAIT","MAIT")
  MAIT.object = signature(MAIT.Object = "MAIT")
  MAIT.object <- methods::new("MAIT")




  if (is.null(data_dir)) {
    stop("No input directory was given")
  }
  if (is.null(preproc_params)) {
    stop("No parameters for preprocessing were included")
  }
  if (is.null(peak_table)) {
    stop("No peak was included")
  }
  print(MAIT::resultsPath(MAIT.object))
  MAIT.object@RawData@parameters@sampleProcessing <- parameters


  class <- list.files(data_dir)
  classNum <- vector(length = length(class))
  fileList <- list.files(path = paste(data_dir, list.files(path = data_dir),
                                      sep = "/"), full.names = TRUE)
  for (i in 1:length(class)) {
    classNum[i] <- length(list.files(paste(data_dir, list.files(data_dir)[i],
                                           sep = "/")))
  }

  classes <- rep(class, classNum)

  if (length(list.files(data_dir)) == 1) {
    warning("Warning: Input data only has one class!")
  }

  if (is.null(project)) {
    warning("Warning: Project name is empty!")
  }

  if (!is.null(project_dir)) {
    resultsPath <- paste(project_dir, paste("results", project, sep = "_"), sep ="/")
    if (any(base::dir.exists(resultsPath))) {
      cat("There are already directories / files in the folder. Not saving new ones.")
      cat("\n")
      }else{
        dir.create(resultsPath)
        peak_table <- base::suppressMessages(as(peak_table, "xcmsSet"))
        peak_table <-  base::suppressMessages(
          base::suppressWarnings(quiet(xcms::fillPeaks(peak_table))
          )
        )
      }
  }else{
        resultsPath <- "results"
        dir.create(resultsPath)
  }

  lcms_write_parameter_table(base::suppressMessages(base::suppressWarnings(
                                                                        quiet(MAIT::parameters(MAIT.object)))),
                          folder = resultsPath)

  peak_table <- base::suppressMessages(as(peak_table, "xcmsSet"))
  peak_table <-  base::suppressMessages(
                      base::suppressWarnings(quiet(xcms::fillPeaks(peak_table))
                                            )
                                   )

  fPeaks <- list(peak_table)
  names(fPeaks) <- "xcmsSet"
  MAIT.object@RawData@data <- fPeaks
  MAIT.object@PhenoData@classes <- class
  MAIT.object@PhenoData@classNum <- classNum
  MAIT.object@PhenoData@resultsPath <- resultsPath
  MAIT <- MAIT.object
}
