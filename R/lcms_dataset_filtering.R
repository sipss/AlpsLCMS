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
#' @return A filtered [lcms_dataset_family] object with the selected polarity
#' @export
#' @examples
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata","dataset_metadata.rds",package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
lcms_filter_polarity <- function(object, polarity.) {
  if (missing(polarity.)) return(object)
  polarity. <- as.numeric(polarity.)
  subset <- MSnbase::polarity(object) %in% polarity.
  object[subset]
  object
}


#' Filter by mass/charge
#'
#' The function filters the dataset by m/z
#' @param dataset An [lcms_dataset_family] object
#' @param mz  The range of masses to filter
#' @family lcms_dataset functions
#' @family lcms_filtering functions
#' @return A filtered [lcms_dataset_family] object with the selecter m/z range
#' @export
#' @examples

#' mz_range = c(300, 400)
#' dataset <- lcms_dataset_load(system.file("extdata","dataset_pos.rds",package = "NIHSlcms"))
#' dataset_pos_mz <-lcms_filter_mz(dataset, mz = mz_range)
#'
#' print(dataset_pos_mz)

lcms_filter_mz <- function(dataset, mz){
      dataset <- MSnbase::filterMz(dataset, mz)

}


#' Filter by retention time
#'
#' The function converts seconds into minutes to cut and keep
#' a range of the retention time in minutes.
#' @param dataset An [lcms_dataset_family] object filtered by retention time
#' @param rt Range of the retention time to keep in minutes
#' @return A filtered dataset with the range of selected rt
#' @export
#' @family lcms_dataset functions
#' @family lcms_filtering functions
#'
#' @export
#' @examples
#'
#' rtime_range = c(5,10)
#' dataset <- lcms_dataset_load(system.file("extdata","dataset_pos.rds",package = "NIHSlcms"))
#' dataset_pos_rt <-lcms_filter_rt_min(dataset, rt = rtime_range)
#'
#' print(dataset_pos_rt)
#'
lcms_filter_rt_min <- function (dataset, rt = c(4, 14)){
  min2sec <- 60
  dataset <- MSnbase::filterRt(dataset, rt = rt * min2sec)
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
#' @param dataset An [lcms_dataset_family] object
#' @param especial_samples A list with the especial samples names.
#' Use `NULL` if there is not any especial sample in the dataset.
#' @return A filtered dataset with the remained samples
#' @export
#' @family lcms_dataset functions
#' @family lcms_filtering functions
#' @export
#'
#' @examples
#' dataset <- lcms_dataset_load(system.file("extdata","dataset_pos_rt.rds",package = "NIHSlcms"))
#' especial_samples <-list(QC = NULL, blank = NULL)
#' datasets_by_class_type <- lcms_filter_sample_type(dataset, especial_samples)
#' dataset_pos_rt_rs <-datasets_by_class_type$regular_samples
#' dataset_pos_rt_qcs <-datasets_by_class_type$QCs
#' print(dataset_pos_rt_rs)
#'
#' print(dataset_pos_rt_qcs)

lcms_filter_sample_type <- function(dataset,  especial_samples){

  QC_label <- especial_samples$QC
  blank_label <- especial_samples$blank

  QC_index <- which(dataset$treatment == QC_label)
  if (length(QC_index) == 0){
    QCs <- NULL
  } else{
    QCs <- MSnbase::filterFile(dataset,file = QC_index)
  }


  blank_index <- which(dataset$treatment == blank_label)
  if (length(blank_index) == 0){
    blanks <- NULL
  } else{
    blanks <- MSnbase::filterFile(dataset,file = blank_index)
  }

  sample_index <- which(!(dataset$treatment %in% c(QC_label, blank_label)))
  if (length(sample_index) == 0){
    regular_samples <- NULL
    stop("Your dataset doesn't have any sample not considered Blank or QC sample")
  } else{
    regular_samples <- MSnbase::filterFile(dataset, file = sample_index)
  }

  datasets_by_class_type <- list(regular_samples = regular_samples,
                                 QCs = QCs, blanks = blanks)
}
