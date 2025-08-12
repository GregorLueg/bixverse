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

# matrix to binary on disk -----------------------------------------------------

set.seed(123)
genes <- sprintf("gene_%i", 1:100)
count_range <- 1:50
no_cells <- 1000L
probs <- 1 / seq_len(length(genes))
probs <- probs / sum(probs)

mirai::daemons(as.integer(parallel::detectCores() / 2))
sparse_data <- mirai::mirai_map(
  1:no_cells,
  \(idx, count_range, genes, probs) {
    set.seed(idx)

    no_expressed <- sample(5:25, 1)

    random_indx_i <- sample(1:length(genes), no_expressed, prob = probs)
    random_counts_i <- sample(count_range, no_expressed, replace = TRUE)

    data.table::data.table(
      i = random_indx_i,
      j = idx,
      x = random_counts_i
    )
  },
  .args = list(count_range, genes, probs)
)[]
mirai::daemons(0)

sparse_dt <- data.table::rbindlist(sparse_data)

csc_matrix <- Matrix::sparseMatrix(
  i = sparse_dt$i,
  j = sparse_dt$j,
  x = sparse_dt$x,
  dims = c(length(genes), no_cells),
  dimnames = list(genes, sprintf("cell_%i", 1:no_cells))
)

## write file ------------------------------------------------------------------

dir <- tempdir()
f_path <- file.path(dir, "test.bin")

single_cell_counts <- SingeCellCountData$new(f_path = f_path)

single_cell_counts$r_csc_mat_to_file(
  no_cells = csc_matrix@Dim[2],
  no_genes = csc_matrix@Dim[1],
  data = as.integer(csc_matrix@x),
  col_ptr = csc_matrix@p,
  row_idx = csc_matrix@i,
  target_size = 1e5
)

expect_true(
  current = "test.bin" %in% list.files(dir),
  info = "test file was written"
)

## get full data back ----------------------------------------------------------

full_data <- single_cell_counts$file_to_r_csc_mat()

expect_equal(
  current = full_data$col_ptr,
  target = csc_matrix@p,
  info = "Recovered column pointers"
)

expect_equal(
  current = full_data$row_idx,
  target = csc_matrix@i,
  info = "Recovered row indices"
)

expect_equal(
  current = full_data$data,
  target = csc_matrix@x,
  info = "Recovered data"
)

## get cells by indices --------------------------------------------------------

indices <- c(1, 2, 10, 15)

r_subsetted_matrix <- csc_matrix[, indices]

rust_subsetted_data <- single_cell_counts$get_cells_by_indices(
  indices = as.integer(indices)
)

expect_equal(
  current = rust_subsetted_data$col_ptr,
  target = r_subsetted_matrix@p,
  info = "Recovered column pointers - subset"
)

expect_equal(
  current = rust_subsetted_data$row_idx,
  target = r_subsetted_matrix@i,
  info = "Recovered row indices - subset"
)

expect_equal(
  current = rust_subsetted_data$data,
  target = r_subsetted_matrix@x,
  info = "Recovered data - subset"
)
