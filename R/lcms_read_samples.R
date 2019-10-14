#' Read LCMS samples on mzXML format
#'
#' @inheritDotParams MSnbase::readMSData
#'
#' @inherit MSnbase::readMSData return
#' @export
#'
#' @examples
#' data(samples_mzxml)
#' dataset <- suppressWarnings(lcms_read_samples(samples_mzxml, mode = "onDisk"))
#' print(dataset)
#'
#'
lcms_read_samples <- function(...){
  dataset <- MSnbase::readMSData(...)
  dataset
}

#
#' Add metadata to MSnExp object
#'
#' @inheritParams Biobase::pData
#' @param metadata A data frame to be merged
#' @param by A column present both in `metadata` and in `Biobase::pData(object)`
#'
#' @return The object with the added metadata
#' @export
#'
#'
lcms_meta_add <- function(object, metadata, by = "sampleNames") {
  phenotype_data <- Biobase::pData(object)
  phenotype_data$sampleNames <- as.character(phenotype_data$sampleNames)
  phenotype_data_extra <- dplyr::left_join(phenotype_data, metadata, by = by)
  Biobase::pData(object) <- phenotype_data_extra
  object
}
