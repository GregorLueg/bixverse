
rextendr::document()

feature_names <- c("A", "B", "C", "D", "E")
total_len <- length(feature_names)
shift <- 1L

feature_a <- purrr::map(1:total_len, \(idx) {
  rep(feature_names[[idx]], total_len - idx + 1 - shift)
})
feature_a <- do.call(c, feature_a)

feature_b <- purrr::map(1:total_len, \(idx) {
  start_point <- idx + shift
  if (start_point <= total_len) feature_names[start_point:total_len] else character(0)
})
feature_b <- do.call(c, feature_b)


rs_diagonal_feature_lists(feature_names, 1L)
