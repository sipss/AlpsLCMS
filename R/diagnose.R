#' Set/Get diagnostic information
#'
#' Set or get diagnostic information
#'
#' @param dataset An [lcms_dataset_family] object
#'
#' @return Diagnostic information, usually a list with data or plots
#' @export
#'
lcms_diagnose <- function(dataset) {
  attr(dataset, "diagnostic")
}


#' @rdname lcms_diagnose
#' @param value The diagnostic we want to set
#' @export
"lcms_diagnose<-" <- function(dataset, value) {
  attr(dataset, "diagnostic") <- value
  dataset
}
