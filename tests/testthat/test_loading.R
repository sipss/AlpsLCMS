context("data-loading")

test_that("loading works", {
  library(faahKO)
  path <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
              recursive = TRUE)[1]
  dataset <- readMSData(path, mode = "onDisk")
  metadata <- data.frame(sampleNames = basename(path),
                         treatment = "ko",
                         stringsAsFactors = FALSE)
  dataset <- lcms_meta_add(dataset, metadata, by = "sampleNames")
  ph <- phData(dataset)[1,2]

  # check result is consistent:
  expect_true(class(dataset@featureData@data[["totIonCurrent"]]) == "numeric")
})

test_that("metadata merge works", {
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
