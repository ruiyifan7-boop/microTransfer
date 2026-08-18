test_that("held-out prediction returns fold-level metrics", {
  dat <- make_test_data()
  fit <- validate_prediction(
    dat$counts,
    dat$metadata,
    "sample_id",
    "study",
    "class",
    features = paste0("f", 9:11)
  )
  expect_s3_class(fit, "mt_prediction_validation")
  expect_equal(nrow(fit$summary), 4)
  expect_true(all(fit$summary$balanced_accuracy >= 0 &
                    fit$summary$balanced_accuracy <= 1))

  null <- prediction_permutation_test(
    dat$counts,
    dat$metadata,
    "sample_id",
    "study",
    "class",
    features = paste0("f", 9:11),
    permutations = 9,
    seed = 2
  )
  expect_equal(nrow(null$null), 9)
  expect_true(null$balanced_accuracy_p > 0 &&
                null$balanced_accuracy_p <= 1)
})
