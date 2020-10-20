context("data-loading")

test_that("readMSData works", {
  library(faahKO)
  path <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
              recursive = TRUE)[1]
  dataset <- MSnbase::readMSData(path, mode = "onDisk")
  dataset
  # check result is consistent:
  expect_true(class(dataset@featureData@data[["totIonCurrent"]]) == "numeric")
  expect_identical(dataset@featureData@data[["retentionTime"]][1], 2501.378)

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

  expect_equal(ph, as.character("ko"))
})
