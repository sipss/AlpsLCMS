context("data-loading")

test_that("readMSData works", {
  library(faahKO)
  path <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
              recursive = TRUE)[1]
  dataset <- readMSData(path, mode = "onDisk")
  # check result is consistent:
  expect_true(is.numeric(dataset@featureData@data[["totIonCurrent"]]))
  expect_true(is.numeric(dataset@featureData@data[["retentionTime"]][1]))

})

context("metadata-merge")

test_that("lcms_meta_add merge works", {
  library(faahKO)
  path <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
              recursive = TRUE)[1]
  dataset <- readMSData(path, mode = "onDisk")
  metadata <- data.frame(sampleNames = basename(path),
                         treatment = "ko",
                         stringsAsFactors = FALSE)
  lcms_meta_add <- function(object, metadata, by = "sampleNames") {
    #making robust the metadata (remove strange characters and separators and numbers as a first characters)
    #Done for treatment, but possibly useful for other variables (check again in the future)
    pattern <- "[\\\"\\s/\\\\,;.:|#@$%&?!*%+-=><^'(){}\\[\\]]+"
    aux_treatment <- stringr::str_replace_all(metadata$treatment,
                                              pattern =pattern,
                                              replacement ="_")
    starts_with_numbers <- stringr::str_detect(aux_treatment,"^[\\d]+")
    for (i in seq_along(aux_treatment)){
      if(starts_with_numbers[i]){
        aux_treatment[i] = paste0("_", aux_treatment[i])
      }
    }
    metadata$treatment <- aux_treatment
    phenotype_data <- Biobase::pData(object)
    phenotype_data$sampleNames <- as.character(phenotype_data$sampleNames)
    phenotype_data_extra <- dplyr::left_join(phenotype_data, metadata, by = by)
    Biobase::pData(object) <- phenotype_data_extra
    object
  }
  dataset <- lcms_meta_add(dataset, metadata, by = "sampleNames")
  ph <- Biobase::pData(dataset)[1,2]
  expect_true(is.integer(dataset@featureData@data[["spectrum"]]))
  expect_true(is.character(class(ph)))
})

context("filter")

test_that("filter works", {
  library(faahKO)
  path <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
              recursive = TRUE)[1]
  dataset <- readMSData(path, mode = "onDisk")
  dataset_shorterRT <- lcms_filter_rt_min(dataset, rt = c(50, 60))
  dataset_shorterMZ <- lcms_filter_mz(dataset_shorterRT, mz = c(200, 500))
  # check result is consistent:
  expect_true(is.numeric(dataset_shorterRT@featureData@data[["totIonCurrent"]]))
  expect_true(is.numeric(dataset_shorterMZ@featureData@data[["totIonCurrent"]]))
})
