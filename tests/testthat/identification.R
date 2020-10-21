context("known_metabolites ")

test_that("known_metabolites works", {
  feature.table <- data.frame(sampleNames = basename(path),
                         treatment = rep("ko",2),
                         stringsAsFactors = FALSE)
  metabolites <- data.frame(sampleNames = basename(path),
                              treatment = rep("ko",2),
                              stringsAsFactors = FALSE)
  # expect_true(is.data.frame(tics))
  # expect_true(is.list(plot))

})
