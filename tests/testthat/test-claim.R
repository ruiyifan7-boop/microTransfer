test_that("claim declaration enforces property-specific outcomes", {
  x <- declare_transferability("membership", "study")
  expect_s3_class(x, "mt_claim")
  expect_error(
    declare_transferability("response", "study"),
    "outcome_col"
  )
})
