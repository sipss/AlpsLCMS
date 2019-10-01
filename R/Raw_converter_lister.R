#' RAW converter
#'
#' The function lists and converts samples in ".raw" format to ".mzXML" format.
#' It is the first step to create a `lcms_dataset`. It requires the previous
#' download and installation of the *RawConverter* application.
#'
#' @param sample_path Directory in which the samples are.
#' @param file_format Format of the LC-MS files (e.g. file_format = "raw").
#' @param rawconverter Directory in where the RawConverter application is located.
#'
#' @return a list of the LC-MS files in a readable format to create the `lcms_dataset`.
#' @export
#'
#' @examples
#' \dontrun{
#' rawconverter_path <-"C:/Users/Biosignal/Desktop/Nestle_2019/LCMS"
#' path <- system.file("extdata", "1.mzXML", package = "NIHSlcms")
#' samples_mzxml <- list_mzxml_samples(path,
#'                                     file_format = "mzXML",
#'                                     rawconverter = rawconverter_path)
#'}

#'
#' \dontrun{
#' path <- system.file("extdata", "1.raw", package = "NIHSlcms")
#' samples_mzxml <- list_mzxml_samples(path,
#'                                     file_format = "raw",
#'                                     rawconverter = rawconverter_path)
#'}

#' \dontrun{
#' samples_mzxml <- list_mzxml_samples(path,
#'                                     file_format = "mzXML",
#'                                    rawconverter = rawconverter_path)
#'}

list_mzxml_samples <- function(sample_path, file_format = "raw", rawconverter){

  if (file_format == "raw") {
    #install_RawConverter(rawconverter)
    samples_raw <- fs::dir_ls(path, glob = "*.raw")
    future::plan(multiprocess)
    lcms_raw_to_mzxml(samples = samples_raw, rawconverter = rawconverter)
    future::plan(sequential)
    sample_names_mzxml <- fs::path_abs(fs::dir_ls(path, glob = "*.mzXML"))
  } else if(file_format == "mzXML") {
    sample_names_mzxml<- fs::path_abs(fs::dir_ls(path, glob = "*.mzXML"))
  } else {
    stop("Not allowed file format. Use only either *.raw or *.mzMXML files")
  }
  sample_names_mzxml
}

