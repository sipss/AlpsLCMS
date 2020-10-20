context("TICs")

test_that("lcms_tics works", {
  library(faahKO)
  path <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
              recursive = TRUE)[1:2]
  dataset <- readMSData(path, mode = "onDisk")
  metadata <- data.frame(sampleNames = basename(path),
                         treatment = rep("ko",2),
                         stringsAsFactors = FALSE)
  dataset <- lcms_meta_add(dataset, metadata, by = "sampleNames")
  dataset <- lcms_filter_polarity(dataset)
  polarity <- 1 # 1 for positive mode, 0 for negative mode
  dataset@featureData@data[["polarity"]] <- rep(polarity, length(dataset@featureData@data[["polarity"]]))
  dataset <- lcms_filter_polarity(dataset, 1)
  tics <- lcms_tics(dataset, treatment = "treatment")
  plot <- lcms_plot_tics(tics,
                 treatment = treatment,
                 plot_type = "spec")
  expect_true(is.data.frame(tics))
  expect_true(is.list(plot))

})
