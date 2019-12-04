#' Read LCMS samples on mzXML format
#'
#' @inheritDotParams MSnbase::readMSData
#'
#' @inherit MSnbase::readMSData return
#' @return An object from the [lcms_dataset_family]
#' @export
#' @examples
#' file_path <- system.file("extdata", package = "NIHSlcms")
#' rawconverter <- NULL
#' file_format <- "mzXML"
#' samples_mzxml <- lcms_list_mzxml_samples(file_path,
#'                                         file_format = file_format,
#'                                         rawconverter = rawconverter)
#' samples_mzxml <- as.character(samples_mzxml)
#' dataset <- suppressWarnings(lcms_read_samples(samples_mzxml, mode = "onDisk"))
#'
#' print(dataset)

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
#' @return An object from the [lcms_dataset_family] with metadata added
#' @export
#' @examples
#' dataset <- lcms_dataset_load(system.file
#'                                   ("extdata","dataset.rds",
#'                                   package = "NIHSlcms"))
#'
#' metadata <- lcms_meta_read(system.file("extdata",
#'                                        "metadata.xlsx",
#'                                        package = "NIHSlcms"))
#'
#' dataset_metadata <- lcms_meta_add(dataset,
#'                                metadata,
#'                                by = "sampleNames")
#' print(dataset_metadata)
lcms_meta_add <- function(object, metadata, by = "sampleNames") {

  #making robust the metadata (remove strange characters and separators and numbers as a first characters)
  #Done for treatment, but possibly useful for other variables (check again in the future)
  pattern <- "[\\\"\\s/\\\\,;.:|#@$%&€¿?¡!*%+-=><^´¨`'(){}\\[\\]]+"
  aux_treatment <- stringr::str_replace_all(metadata$treatment,
                                   pattern =pattern,
                                   replacement ="_")
  starts_with_numbers <- stringr::str_detect(aux_treatment,"^[\\d]+")
  for (i in seq_along(aux_treatment)){
    if(starts_with_numbers[i]){
      aux_treatment[i] = paste0("_", aux_treatment[i])
    }
  }
  metadata$treatment <- aux_treatment
  phenotype_data <- Biobase::pData(object)
  phenotype_data$sampleNames <- as.character(phenotype_data$sampleNames)
  phenotype_data_extra <- dplyr::left_join(phenotype_data, metadata, by = by)
  Biobase::pData(object) <- phenotype_data_extra
  object
}


