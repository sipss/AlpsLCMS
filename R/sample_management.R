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

#' Function to load a lcms_dataset object
#'
#' @param file_name The file name to load
#' @return An object from the [lcms_dataset_family]
#' @export
#' @family lcms_dataset functions
#' @family lcms_dataset_peak_table functions
#' @family import/export functions
#'
#' @examples
#' file_name <- system.file("extdata", "dataset.rds", package = "NIHSlcms")
#' dataset <- lcms_dataset_load(file_name)
#' print(dataset)

lcms_dataset_load <- function(file_name) {
  dataset <-base::readRDS(file_name)
  dataset
}


#' Function to save a lcms_dataset object
#'
#' @param dataset An object from the [lcms_dataset_family]
#' @param file_name The file name to save to
#' @param ... Additional arguments passed to [saveRDS].
#' @export
#' @examples
#' dataset <- lcms_dataset_load(system.file("extdata", "dataset.rds", package = "NIHSlcms"))
#' file_name <- "dataset.rds"
#' lcms_dataset_save(dataset, file_name)
#' print(dataset)
#'
lcms_dataset_save <- function(dataset, file_name, ...) {
  lcms_diagnose(dataset) <- NULL
  saveRDS(dataset, file_name)
}

#NULL

#' Export Metadata to an Excel file
#'
#' @param dataset An [lcms_dataset_family] object
#' @param xlsx_file "The .xlsx excel file"
#' @param groups A character vector. Use `"external"` for the external metadata or
#'  the default for a more generic solution
#' @return The Excel file name
#' @export
#' @family metadata functions
#' @family lcms_dataset functions
#' @family lcms_dataset_peak_table functions
#' @family import/export functions
#'
#' @examples
#' dataset_metadata <- lcms_dataset_load(system.file("extdata",
#'                                              "dataset_metadata.rds",
#'                                               package = "NIHSlcms"))
#' xlsx_file <- paste0(system.file("extdata",package = "NIHSlcms"),
#'                     "/", "exported_metadata.xlsx")
#'
#' lcms_meta_export(dataset_metadata, xlsx_file)
#' print(dataset_metadata)
#'
lcms_meta_export <- function(dataset,
                             xlsx_file) {
  groups_present <- phData(dataset)
  writexl::write_xlsx(x = NIHSlcms::phData(dataset), path = xlsx_file)
}

#' Reads metadata from an Excel file
#'
#' @param xlsx_file xlsx_file "The .xlsx excel file" with metadata
#' @return A dataframe with the metadata
#' @export
#' @family metadata functions
#' @family lcms_dataset functions
#' @family lcms_dataset_peak_table functions
#' @family import/export functions
#'
#' @examples
#' metadata <- lcms_meta_read(system.file("extdata",
#'                                              "metadata.xlsx",
#'                                               package = "NIHSlcms"))
#' print(metadata[1:6])
#'
lcms_meta_read <- function(xlsx_file) {
  meta <- readxl::read_excel(xlsx_file)
  meta
}

#NULL

#' Generic function to access the phenotypic data (Biobase)
#'
#' It uses [Biobase::pData] to access phenotypic data such as
#' metadata, covariates, etc...
#'
#' @param object a `lcms_dataset`
#' @inheritParams Biobase::pData
#' @return phenoData returns an object containing information on
#' both variable values and variable meta-data. varLabels returns
#' a character vector of measured variables. pData returns a data
#' frame with samples as rows, variables as columns. varMetadata
#' returns a data frame with variable names as rows, description
#' tags (e.g., unit of measurement) as columns.
#' @export
#'
#' @examples
#' dataset <- lcms_dataset_load(system.file("extdata",
#'                                          "dataset.rds",
#'                                           package = "NIHSlcms"))
#' xlsx_file <- system.file("extdata", "exported_metadata.xlsx", package = "NIHSlcms")
#' metadata<- readxl::read_excel(xlsx_file)
#' dataset_metadata <- lcms_meta_add(dataset, metadata, by = "sampleNames")
#'
#' head(phData(dataset_metadata))[1:4]

phData <- function (object) {
  Biobase::pData(object)
}

