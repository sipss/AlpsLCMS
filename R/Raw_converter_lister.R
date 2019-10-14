#' RAW converter
#'
#' The function lists and converts samples in ".raw" format to ".mzXML" format.
#' It is the first step to create a `lcms_dataset`. It requires the previous
#' download and installation of the *RawConverter* application.
#'
#' @param sample_path Directory in which the samples are.
#' @param file_format Format of the LC-MS files (e.g. file_format = "raw").
#' @param rawconverter_path Directory in where the RawConverter application is located.
#'
#' @return A list of the LC-MS files in a readable format to create the `lcms_dataset`.
#' @export
#'
#' @examples
#' sample_path <- system.file("extdata", package = "NIHSlcms")
#' samples_mzxml <- lcms_list_mzxml_samples(sample_path,
#'                                     file_format = "mzXML")
#'
#' \dontrun{
#' samples_mzxml <- lcms_list_mzxml_samples(sample_path,
#'                                     file_format = "raw",
#'                                     rawconverter_path = rawconverter_path)
#'}

lcms_list_mzxml_samples <- function(sample_path,
                                    file_format = "mzXML",
                                    rawconverter_path = NULL){

  if (file_format == "raw") {
    #install_RawConverter(rawconverter)
    samples_raw <- fs::dir_ls(sample_path, glob = "*.raw")
    future::plan(multiprocess)
    lcms_raw_to_mzxml(samples = samples_raw, rawconverter = rawconverter_path)
    future::plan(sequential)
    sample_names_mzxml <- fs::path_abs(fs::dir_ls(sample_path, glob = "*.mzXML"))
  } else if(file_format == "mzXML") {
    sample_names_mzxml<- fs::path_abs(fs::dir_ls(sample_path, glob = "*.mzXML"))
  } else {
    stop("Not allowed file format. Use only either *.raw or *.mzMXML files")
  }
  sample_names_mzxml
}

