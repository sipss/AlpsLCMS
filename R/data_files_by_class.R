#' Title
#'
#' @param x
#' @param lcms_dataset An [lcms_dataset_family] object
#' @return
#' @export
#'
#' @examples
rearrange_datafiles_by_class = function(lcms_dataset, dataDir) {
  no_blank_files <- pData(lcms_dataset)$sampleNames
  no_blank_files_treatment <- pData(lcms_dataset)$treatment

  filetreat_info <- data.frame(Filename = no_blank_files,
                               Treatment= no_blank_files_treatment,
                               stringsAsFactors = FALSE)
  filetreat_info <- filetreat_info %>% group_by(Treatment) %>% arrange(Treatment)
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
