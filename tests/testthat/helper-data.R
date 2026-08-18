make_test_data <- function() {
  set.seed(42)
  metadata <- data.frame(
    sample_id = paste0("s", 1:48),
    study = rep(paste0("C", 1:4), each = 12),
    x = rep(seq(-1, 1, length.out = 12), 4),
    class = rep(rep(c("A", "B"), each = 6), 4),
    stringsAsFactors = FALSE
  )
  counts <- matrix(
    rpois(48 * 20, 0.2),
    nrow = 48,
    dimnames = list(metadata$sample_id, paste0("f", 1:20))
  )
  counts[, 1:3] <- counts[, 1:3] + 5
  counts[metadata$study == "C1", 4:8] <-
    counts[metadata$study == "C1", 4:8] + 6
  counts[metadata$class == "B", 9:11] <-
    counts[metadata$class == "B", 9:11] + 4
  list(counts = counts, metadata = metadata)
}
