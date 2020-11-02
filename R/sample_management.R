#' Read mzXML samples
#'
#' Read LCMS samples on mzXML format.
#'
#' @inheritDotParams MSnbase::readMSData
#' @inherit MSnbase::readMSData return
#' @return A lcms_dataset.
#' @family import/export functions
#' @export
#' @examples
#' \dontrun{
#' file_path <- system.file("extdata", package = "AlpsLCMS")
#' rawconverter <- NULL
#' file_format <- "mzXML"
#' samples_mzxml <- lcms_list_mzxml_samples(file_path,
#'                                         file_format = file_format,
#'                                         rawconverter = rawconverter)
#' samples_mzxml <- as.character(samples_mzxml)
#' dataset <- suppressWarnings(lcms_read_samples(samples_mzxml, mode = "onDisk"))
#'
#' print(dataset)
#' }

lcms_read_samples <- function(...){
  dataset <- MSnbase::readMSData(...)
  dataset
}

#' Add metadata
#'
#' Add metadata to MSnExp object.
#'
#' @param object A lcms_dataset.
#' @param metadata A data frame to be merged.
#' @param by A column present both in `metadata` and in `Biobase::pData(object)`.
#' @return A lcms_dataset with the added metadata.
#' @family metadata functions
#' @family dataset functions
#' @export
#' @examples
#' \dontrun{
#' dataset <- lcms_dataset_load(system.file
#'                                   ("extdata","dataset.rds",
#'                                   package = "AlpsLCMS"))
#'
#' metadata <- lcms_meta_read(system.file("extdata",
#'                                        "metadata.xlsx",
#'                                        package = "AlpsLCMS"))
#'
#' dataset_metadata <- lcms_meta_add(dataset,
#'                                metadata,
#'                                by = "sampleNames")
#' print(dataset_metadata)
#' }
lcms_meta_add <- function(object, metadata, by = "sampleNames") {

  #making robust the metadata (remove strange characters and separators and numbers as a first characters)
  #Done for treatment, but possibly useful for other variables (check again in the future)
  pattern <- "[\\\"\\s/\\\\,;.:|#@$%&?!*%+-=><^'(){}\\[\\]]+"
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

#' Load lcms_datasets
#'
#' Function to load a lcms_dataset object.
#'
#' @param file_name The file name to load.
#' @return A lcms_dataset.
#' @family dataset functions
#' @family dataset_peak_table functions
#' @family import/export functions
#' @export
#' @examples
#' \dontrun{
#' file_name <- system.file("extdata", "dataset.rds", package = "AlpsLCMS")
#' dataset <- lcms_dataset_load(file_name)
#' print(dataset)
#' }

lcms_dataset_load <- function(file_name) {
  dataset <-base::readRDS(file_name)
  dataset
}


#' Save lcms_datasets
#'
#' Function to save a lcms_dataset object.
#'
#' @param dataset A lcms_dataset.
#' @param file_name The file name to save to.
#' @param ... Additional arguments passed to [saveRDS].
#' @family dataset functions
#' @family import/export functions
#' @export
#' @examples
#' \dontrun{
#' dataset <- lcms_dataset_load(system.file("extdata", "dataset.rds", package = "AlpsLCMS"))
#' file_name <- "dataset.rds"
#' lcms_dataset_save(dataset, file_name)
#' print(dataset)
#' }
lcms_dataset_save <- function(dataset, file_name, ...) {
  lcms_diagnose(dataset) <- NULL
  saveRDS(dataset, file_name)
}


#' Export metadata
#'
#' Export Metadata to an Excel file.
#'
#' @param dataset A lcms_dataset.
#' @param xlsx_file The .xlsx excel file.
#' @return The Excel file name.
#' @family metadata functions
#' @family dataset functions
#' @family import/export functions
#' @export
#' @examples
#' \dontrun{
#' dataset_metadata <- lcms_dataset_load(system.file("extdata",
#'                                              "dataset_metadata.rds",
#'                                               package = "AlpsLCMS"))
#' xlsx_file <- paste0(system.file("extdata",package = "AlpsLCMS"),
#'                     "/", "exported_metadata.xlsx")
#'
#' lcms_meta_export(dataset_metadata, xlsx_file)
#' print(dataset_metadata)
#' }
lcms_meta_export <- function(dataset,
                             xlsx_file) {
  groups_present <- phData(dataset)
  writexl::write_xlsx(x = AlpsLCMS::phData(dataset), path = xlsx_file)
}

#' Read metadata
#'
#' Reads metadata from an Excel file.
#'
#' @param xlsx_file xlsx_file "The .xlsx excel file" with metadata.
#' @return A dataframe with the metadata.
#' @family metadata functions
#' @family import/export functions
#' @export
#' @examples
#' \dontrun{
#' metadata <- lcms_meta_read(system.file("extdata",
#'                                              "metadata.xlsx",
#'                                               package = "AlpsLCMS"))
#' print(metadata[1:6])
#' }
lcms_meta_read <- function(xlsx_file) {
  meta <- readxl::read_excel(xlsx_file)
  meta
}


#' Generic function to access the phenotypic data (Biobase)
#'
#' It uses [Biobase::pData] to access phenotypic data such as
#' metadata, covariates, etc...
#' @inheritParams Biobase::pData
#' @return phenoData returns an object containing information on
#' both variable values and variable meta-data. varLabels returns
#' a character vector of measured variables. pData returns a data
#' frame with samples as rows, variables as columns. varMetadata
#' returns a data frame with variable names as rows, description
#' tags (e.g., unit of measurement) as columns.
#' @family metadata functions
#' @family dataset functions
#' @family import/export functions
#' @export
#' @examples
#' \dontrun{
#' dataset <- lcms_dataset_load(system.file("extdata",
#'                                          "dataset.rds",
#'                                           package = "AlpsLCMS"))
#' xlsx_file <- system.file("extdata", "exported_metadata.xlsx", package = "AlpsLCMS")
#' metadata<- readxl::read_excel(xlsx_file)
#' dataset_metadata <- lcms_meta_add(dataset, metadata, by = "sampleNames")
#'
#' head(phData(dataset_metadata))[1:4]
#' }
phData <- function (object) {
  Biobase::pData(object)
}

#' is.negative
#'
#' Function to confirm if the dataset was acquired in negative mode.
#'
#' @param dataset
#'
#' @return
#' @export
is.negative <- function (dataset){
  negative_polarity = c(0,0,0,0)
  actual_polarity = dataset@featureData@data[["polarity"]][1:4]
  message("dataset polarity is negative: ", all.equal(actual_polarity, negative_polarity))
}

#' is.positive
#'
#' Function to confirm if the dataset was acquired in positive mode.
#'
#' @param dataset
#'
#' @return
#' @export
is.positive <- function (dataset){
  positive_polarity = c(1,1,1,1)
  actual_polarity = dataset@featureData@data[["polarity"]][1:4]
  message("dataset polarity is positive: ", all.equal(actual_polarity, positive_polarity))
}
