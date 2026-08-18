test_that("matched null sampling preserves size, uniqueness, and strata", {
  universe <- paste0("f", seq_len(100))
  strata <- setNames(rep(c("A", "B", "C", "D"), each = 25), universe)
  core <- c(universe[1:5], universe[26:30], universe[51:55])

  set.seed(7)
  sampled <- .sample_matched_set(core, universe, strata)

  expect_length(sampled, length(core))
  expect_length(unique(sampled), length(core))
  expect_false(any(sampled %in% core))
  expect_equal(
    sort(table(strata[sampled])),
    sort(table(strata[core]))
  )
})
