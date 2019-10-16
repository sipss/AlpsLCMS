#' Filter an experiment by its polarity
#'
#' This function was merged in
#' the MSnbase package on 2018-12-10 as proposed at:
#'
#' https://github.com/lgatto/MSnbase/issues/388
#'
#' Once MSnbase is released (Bioconductor 3.9, maybe once R 3.6 is out),
#' remove this function.
#'
#'
#' @param object An MSnExp object
#' @param polarity. The polarity to keep
#' @family lcms_dataset functions
#' @family lcms_filtering functions
#' @export
#' @examples
#'
#' lcms_dataset_2_polarities <- lcms_dataset_load(system.file("extdata","lcms_dataset_metadata.rds",package = "NIHSlcms"))
#' lcms_dataset_pos <- lcms_filterPolarity(lcms_dataset_2_polarities, polarity. = 1)
#'
lcms_filterPolarity <- function(object, polarity.) {
  if (missing(polarity.)) return(object)
  polarity. <- as.numeric(polarity.)
  subset <- MSnbase::polarity(object) %in% polarity.
  object[subset]
  object
}

#' Filter by retention time
#'
#' The function converts seconds into minutes to cut and keep
#' a range of the retention time in minutes.
#' @param lcms_dataset An [lcms_dataset_family] object
#' @param rt Range of the retention time to keep in minutes
#' @return A filtered dataset with the range of selected rt
#' @export
#' @family lcms_dataset functions
#' @family lcms_filtering functions
#'
#' @export
#' @examples
#'
#' rtime_range = c(4,14)
#' lcms_dataset <- lcms_dataset_load(system.file("extdata","lcms_dataset_metadata.rds",package = "NIHSlcms"))
#' lcms_dataset_rt <-lcms_filterRTmin(lcms_dataset, rt = rtime_range)
#'
lcms_filterRTmin <- function (lcms_dataset, rt = c(4, 14)){
  min2sec <- 60
  lcms_dataset <- MSnbase::filterRt(lcms_dataset, rt = rt * min2sec)
}

#' Filter by sample treatment
#'
#' In a dataset, there are different types of samples for
#' checking purposes. For instance, quality control (QC) samples
#' may be pools of all samples, and blank samples may be composed
#' by the solvent used. Therefore, these samples should be removed
#' before alignment and processing from the main dataset. Use `NULL`
#' if there is not any especial sample in the dataset.
#'
#' The function can distinguish and filter different samples types:
#' * Regular samples
#' * Blank samples
#' * Quality Control samples
#' @param lcms_dataset An [lcms_dataset_family] object
#' @param especial_samples A list with the especial samples names.
#' Use `NULL` if there is not any especial sample in the dataset.
#' @return A filtered dataset with the remained samples
#' @export
#' @family lcms_dataset functions
#' @family lcms_filtering functions
#' @export
#'
#' @examples

#' lcms_dataset <- lcms_dataset_load(system.file("extdata","lcms_dataset_metadata.rds",package = "NIHSlcms"))
#' especial_samples <-list(QC = "resveratrol", blank = "blank")
#' datasets_by_class_type <- lcms_filterSampleType(lcms_dataset, especial_samples)
#' lcms_dataset_regular_samples <-datasets_by_class_type$regular_samples
#' lcms_dataset_qcs <-datasets_by_class_type$QCs

lcms_filterSampleType <- function(lcms_dataset,  especial_samples){

  QC_label <- especial_samples$QC
  blank_label <- especial_samples$blank

  QC_index <- which(lcms_dataset$treatment == QC_label)
  if (length(QC_index) == 0){
    QCs <- NULL
  } else{
    QCs <- MSnbase::filterFile(lcms_dataset,file = QC_index)
  }


  blank_index <- which(lcms_dataset$treatment == blank_label)
  if (length(blank_index) == 0){
    blanks <- NULL
  } else{
    blanks <- MSnbase::filterFile(lcms_dataset,file = blank_index)
  }

  sample_index <- which(!(lcms_dataset$treatment %in% c(QC_label, blank_label)))
  if (length(sample_index) == 0){
    regular_samples <- NULL
    stop("Your dataset doesn't have any sample not considered Blank or QC sample")
  } else{
    regular_samples <- MSnbase::filterFile(lcms_dataset, file = sample_index)
  }

  datasets_by_class_type <- list(regular_samples = regular_samples,
                                 QCs = QCs, blanks = blanks)
}
