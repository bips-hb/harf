library(testthat)
library(harf)

test_that("kich dataset loads correctly", {
  data(kich)
  expect_true(is.data.frame(kich))
  expect_true(all(c("tumor_stage", "gender") %in% colnames(kich)))
})

test_that("Adversarial game for clustering works crrectly", {
  iris_data <- datasets::iris
    harf_cl_obj <- h_arf(omx_data = iris_data[ , -5],
                      cli_lab_data = iris_data[ , 5, drop = FALSE],
                      chunck_size = 2,
                      target = NULL,
                      parallel = FALSE,
                      verbose = FALSE)
  expect_true(inherits(harf_cl_obj, "harf"))
  # Forge
  synth_iris <- h_forge(
    harf_obj = harf_cl_obj,
    n_synth = 150,
    evidence = NULL,
    verbose = FALSE,
    parallel = FALSE
  )
  expect_true(is.data.frame(synth_iris))
})


test_that("Adversarial game for prediction works correctly", {
  iris_data <- datasets::iris
    harf_prd_obj <- h_arf(omx_data = iris_data[ , -5],
                      cli_lab_data = iris_data[ , 5, drop = FALSE],
                      chunck_size = 2,
                      target = "Species",
                      parallel = FALSE,
                      verbose = FALSE)
  expect_true(inherits(harf_prd_obj, "harf"))
  # Forge
  synth_iris <- h_forge(
    harf_obj = harf_prd_obj,
    n_synth = 150,
    evidence = NULL,
    parallel = FALSE
  )
  expect_true(is.data.frame(synth_iris))
})
