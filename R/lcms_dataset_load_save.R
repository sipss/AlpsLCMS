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
lcms_dataset_save <- function(lcms_dataset, file_name, ...) {
  lcms_diagnose(lcms_dataset) <- NULL
  saveRDS(lcms_dataset, file_name)
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
#' \dontrun{
#' dataset_metadata <- lcms_dataset_load(system.file("extdata",
#'                                              "dataset_metadata.rds",
#'                                               package = "NIHSlcms"))
#' xlsx_file <- "exported_metadata.xlsx"
#' lcms_meta_export(dataset_metadata, xlsx_file)
#' print(dataset_metadata)
#'
#' }
lcms_meta_export <- function(dataset,
                             xlsx_file) {
  groups_present <- phData(lcms_dataset)
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
