test_that("LOSO membership identifies stable features", {
  dat <- make_test_data()
  fit <- validate_membership(
    dat$counts,
    dat$metadata,
    "sample_id",
    "study",
    prevalence = 0.5,
    n_null = 19,
    seed = 1
  )
  expect_s3_class(fit, "mt_membership_validation")
  expect_true(all(c("f1", "f2", "f3") %in% fit$stable_features))
  expect_equal(nrow(fit$summary), 4)
  contrast <- pooled_membership_contrast(
    dat$counts, dat$metadata, "sample_id", "study", 0.5
  )
  expect_equal(nrow(contrast), 1)
  expect_gte(contrast$pooled_candidates,
             contrast$study_aware_candidates)
})
