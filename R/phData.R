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
#' \dontrun{
#' metadata<- readxl::read_excel("D:/myFile.xlsx")
#' dataset <- lcms_meta_add(dataset, metadata, by = "sampleNames")
#' phData(dataset)
#' }

phData <- function (object) {
  Biobase::pData(object)
}
