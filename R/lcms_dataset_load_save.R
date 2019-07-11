#' Functions to load and save lcms_dataset objects
#'
#' @name load_and_save_functions
#' @param file_name The file name to load or save to
#' @param lcms_dataset An object from the [lcms_dataset_family]
#' @param ... Additional arguments passed to [saveRDS].
#' @family lcms_dataset functions
#' @family lcms_dataset_peak_table functions
#' @family import/export functions
NULL

#' @rdname load_and_save_functions
#' @export
lcms_dataset_load <- function(file_name) {
  return(readRDS(file_name))
}

#' @rdname load_and_save_functions
#' @export
lcms_dataset_save <- function(lcms_dataset, file_name, ...) {
  nmr_diagnose(lcms_dataset) <- NULL
  saveRDS(lcms_dataset, file_name)
  return(lcms_dataset)
}
