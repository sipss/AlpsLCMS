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
#' file_name <- system.file("extdata", "lcms_dataset.rds", package = "NIHSlcms")
#' lcms_dataset <- lcms_dataset_load(file_name)
#' print(lcms_dataset)
#'
lcms_dataset_load <- function(file_name) {
  dataset <-base::readRDS(file_name)
  dataset
}


#' Function to save a lcms_dataset object
#'
#' @param lcms_dataset An object from the [lcms_dataset_family]
#' @param file_name The file name to save to
#' @param ... Additional arguments passed to [saveRDS].
#' @export
#' @examples
#'
#' lcms_dataset <- lcms_dataset_load(system.file("extdata", "lcms_dataset.rds", package = "NIHSlcms"))
#' file_name <- "lcms_dataset.rds"
#' lcms_dataset_save(lcms_dataset, file_name)
#' print(lcms_dataset)
#'
lcms_dataset_save <- function(lcms_dataset, file_name, ...) {
  lcms_diagnose(lcms_dataset) <- NULL
  saveRDS(lcms_dataset, file_name)
}

NULL

#' Export Metadata to an Excel file
#'
#' @param lcms_dataset An [lcms_dataset_family] object
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
#' \dontrun{
#' lcms_dataset_metadata <- lcms_dataset_load(system.file("extdata",
#'                                              "lcms_dataset_metadata.rds",
#'                                               package = "NIHSlcms"))
#' xlsx_file <- "exported_metadata.xlsx"
#' lcms_meta_export(lcms_dataset_metadata, xlsx_file)
#' print(lcms_dataset_metadata)
#' }
#'
lcms_meta_export <- function(lcms_dataset,
                             xlsx_file) {
  groups_present <- phData(lcms_dataset)
  writexl::write_xlsx(x = NIHSlcms::phData(lcms_dataset), path = xlsx_file)
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
#'
#' metadata <- lcms_meta_read(system.file("extdata",
#'                                              "metadata.xlsx",
#'                                               package = "NIHSlcms"))
#' print(metadata)
#'
  lcms_meta_read <- function(xlsx_file) {
  meta <- readxl::read_excel(xlsx_file)
  meta
}

NULL
