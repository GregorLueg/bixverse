# libraries --------------------------------------------------------------------

# data -------------------------------------------------------------------------

data <- c(0.3, 0.0, 1.0, 0, 0.2, -0.5)

data_dense <- rs_upper_triangle_to_dense(data, 1, 4)

expected_sparse_matrix <- as(
  Matrix::Matrix(data_dense, sparse = TRUE),
  "generalMatrix"
)

## rust tests ------------------------------------------------------------------

rs_data <- rs_upper_triangle_to_sparse(data, 1, 4)

expect_equal(
  current = rs_data$row_indices,
  target = expected_sparse_matrix@i,
  info = "Rust sparse implementation - row indices"
)

expect_equal(
  current = rs_data$col_ptr,
  target = expected_sparse_matrix@p,
  info = "Rust sparse implementation - col pointers"
)

expect_equal(
  current = unlist(rs_data$data),
  target = expected_sparse_matrix@x,
  info = "Rust sparse implementation - data"
)

## r tests ---------------------------------------------------------------------

sparse_matrix <- upper_triangle_to_sparse(data, 1L, 4L)

expect_equal(
  current = sparse_matrix,
  target = sparse_matrix,
  info = "Upper triangle to sparse via Rust"
)
