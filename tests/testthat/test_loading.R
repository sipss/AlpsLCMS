context("data-loading")

test_that("readMSData works", {
  library(faahKO)
  path <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
              recursive = TRUE)[1]
  dataset <- readMSData(path, mode = "onDisk")

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
  dataset <- lcms_meta_add(dataset, metadata, by = "sampleNames")
  ph <- phData(dataset)[1,2]

  expect_equal(ph, as.character("ko"))
})
