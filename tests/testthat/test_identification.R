context("known_metabolites ")

test_that("known_metabolites works", {
  feature.table <- data.frame(feature = c("one", "two"),
                              mz = as.numeric(c(132.1021, 118.0862)),
                              rt = as.numeric(c(5.55 , 7.77)),
                              stringsAsFactors = FALSE)
  # metabolites <- data.frame(metabolite = c("metab01", "metabo2", "metabo3"),
  #                           mz = as.numeric(c(123.333, 125.555, 332.222)),
  #                           rt = as.numeric(c(5.55, 7.77, 9)),
  #                           stringsAsFactors = FALSE)
  met <- assignation_pos_HMDB(feature.table)

  expect_true(is.character(class(met[length(met), 1])))

})
