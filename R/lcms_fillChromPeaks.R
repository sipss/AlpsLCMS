#' Filling missing values in a peak tabe
#'
#' In the imputation stage, we integrate the areas of the missing peaks of the peak table
#' that were not detected in the previous steps of the signal preprocessing workflow.
#'
#' @param peak_table A table of peaks with (possibly) missing values.
#' @return A peak table where the missing peaks have been filled
#' @export
#' @examples
#' \dontrun{
#' file_name <-  system.file("extdata", "peak_table.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' peak_table_imputed <- lcms_fillChromPeaks(peak_table)
#'
#' print(peak_table_imputed)
#' }
lcms_fillChromPeaks <- function(peak_table){
  quiet <- function(x) {
    base::sink(base::tempfile())
    base::on.exit(base::sink())
    base::invisible(base::force(x))
  }

  cat("\n","Imputing missing peaks of the peak table.","\n")

  peak_table_imputed <-  base::suppressWarnings(
                            base::suppressMessages(quiet(xcms::fillChromPeaks(peak_table))
                                                  )
                                               )
  return(peak_table_imputed)
}
