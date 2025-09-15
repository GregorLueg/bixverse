# simple sparse tests ----------------------------------------------------------

## data ------------------------------------------------------------------------

data <- c(0.3, 0.0, 1.0, 0, 0.2, -0.5)

data_dense <- rs_upper_triangle_to_dense(data, 1, 4)

expected_sparse_matrix <- as(
  Matrix::Matrix(data_dense, sparse = TRUE),
  "generalMatrix"
)

### rust tests -----------------------------------------------------------------

rs_data <- rs_upper_triangle_to_sparse(data, 1, 4)

sparse_matrix <- sparse_list_to_mat(rs_data)

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
  current = rs_data$data,
  target = expected_sparse_matrix@x,
  info = "Rust sparse implementation - data"
)

expect_true(
  current = rs_data$ncol == 4,
  info = "Rust sparse implementation - ncols"
)

expect_true(
  current = rs_data$nrow == 4,
  info = "Rust sparse implementation - nrow"
)

expect_equal(
  current = sparse_matrix,
  target = expected_sparse_matrix,
  info = "Rust sparse list to sparse mat"
)

### r tests --------------------------------------------------------------------

sparse_matrix <- upper_triangle_to_sparse(data, 1L, 4L)

expect_equal(
  current = sparse_matrix,
  target = expected_sparse_matrix,
  info = "Upper triangle to sparse via Rust"
)

### class test -----------------------------------------------------------------

rownames(expected_sparse_matrix) <- colnames(expected_sparse_matrix) <- c(
  "a",
  "b",
  "c",
  "d"
)

upper_triangle_repr <- bixverse:::upper_triangular_sym_mat$new(
  values = data,
  features = c("a", "b", "c", "d"),
  shift = 1L
)

sparse_matrix_class <- upper_triangle_repr$get_sparse_matrix(.verbose = FALSE)

expect_equal(
  current = sparse_matrix_class,
  target = expected_sparse_matrix,
  info = "Upper triangle to sparse via Rust"
)
