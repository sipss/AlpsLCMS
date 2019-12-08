#' Set/Get diagnostic information
#'
#' Set or get diagnostic information
#'
#' @param lcms_dataset An [lcms_dataset_family] object
#'
#' @return Diagnostic information, usually a list with data or plots
#' @export
#'
lcms_diagnose <- function(lcms_dataset) {
  attr(lcms_dataset, "diagnostic")
}


#' @rdname lcms_diagnose
#' @param value The diagnostic we want to set
#' @export
"lcms_diagnose<-" <- function(lcms_dataset, value) {
  attr(lcms_dataset, "diagnostic") <- value
  lcms_dataset
}
