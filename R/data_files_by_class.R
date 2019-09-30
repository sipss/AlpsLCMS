#' Rearrange datafiles by class
#'
#' Although Data Preprocessing can be performed using any
#' filepath with the `XCMS` Package, it is convenient to
#' rearrange data files by sample class (that is, all the
#' samples belonging to the same sample class are included
#' in the same folder). We do this because, after the
#' Preprocessing stage, the `MAIT` Package takes care
#' of Data Annotation and Metabolite Identification. `MAIT`
#' needs a specific directory structure for data managing.
#' Otherwise, it can't work properly.
#'
#' @param lcms_dataset An [lcms_dataset_family] object
#' @param dataDir a directory in where LC-MS files are
#' going to be saved
#' @return LC-MS datafiles sorted by class treatment
#' @export
#'
#' @examples
#' rearrange_datafiles_by_class(lcms_dataset = lcms_dataset,
#'                              dataDir = path)
rearrange_datafiles_by_class = function(lcms_dataset, dataDir) {
  no_blank_files <- Biobase::pData(lcms_dataset)$sampleNames
  no_blank_files_treatment <- Biobase::pData(lcms_dataset)$treatment

  filetreat_info <- data.frame(Filename = no_blank_files,
                               Treatment= no_blank_files_treatment,
                               stringsAsFactors = FALSE)
  filetreat_info <- filetreat_info %>% dplyr::group_by(Treatment) %>% dplyr::arrange(Treatment)
  filetreat_info <- split(filetreat_info, filetreat_info$Treatment)

  drop_treatment <-function(x) {
    x$Treatment <- NULL
    colnames(x) <- "FileName"
    x
  }

  filetreat_info <- lapply(filetreat_info, FUN = drop_treatment)

  for (i in seq_along(filetreat_info)){
  filer <- filetreat_info[[i]][["FileName"]]
  foldname <-names(filetreat_info)[i]
  treatDir <- paste0(dataDir, "/", foldname, "/")
  dir.create(treatDir)
  data_subset <- lcms_dataset %>% MSnbase::filterFile(file = filer)
  Biobase::fData(data_subset)$centroided <- TRUE
  Biobase::fData(data_subset)$peaksCount <- Biobase::fData(data_subset)$originalPeaksCount
  mzR::writeMSData(data_subset,
                   file = unlist(lapply(treatDir, FUN = paste0, filer)),
                   outformat = c("mzxml"), copy = FALSE)
  }
}
