# sparse data format -----------------------------------------------------------

library(magrittr)

## csc (cell-centric) ----------------------------------------------------------

set.seed(123)
genes <- sprintf("gene_%i", 1:20000)
count_range <- 1:20
no_cells <- 10000L
probs <- 1 / seq_len(length(genes))
probs <- probs / sum(probs)

# More efficient approach: build sparse matrix directly
mirai::daemons(5L)

# Generate all random data in parallel chunks
sparse_data <- mirai::mirai_map(
  1:no_cells,
  \(idx, count_range, genes, probs) {
    set.seed(idx)

    random_indx_i <- sample(1:length(genes), 750, prob = probs)
    random_counts_i <- sample(count_range, 750, replace = TRUE)

    data.table::data.table(
      i = random_indx_i,
      j = idx,
      x = random_counts_i
    )
  },
  .args = list(count_range, genes, probs)
)[]

sparse_dt <- data.table::rbindlist(sparse_data)

csc_matrix <- Matrix::sparseMatrix(
  i = sparse_dt$i,
  j = sparse_dt$j,
  x = sparse_dt$x,
  dims = c(length(genes), no_cells),
  dimnames = list(genes, sprintf("cell_%i", 1:no_cells))
)

mirai::daemons(0L) # Clean up daemons

rextendr::document()

csc_matrix@i
csc_matrix@p
csc_matrix@x
csc_matrix@Dim

length(csc_matrix@p)

dir <- tempdir()
f_path <- file.path(dir, "test.bin")

tictoc::tic()
rs_csc_to_binary_f(
  f_path = f_path,
  no_cells = csc_matrix@Dim[2],
  no_genes = csc_matrix@Dim[1],
  data = as.integer(csc_matrix@x),
  col_ptr = csc_matrix@p,
  row_idx = csc_matrix@i
)
tictoc::toc()

list.files(dir)

tictoc::tic()
read_file <- rs_binary_f_to_csc(f_path)
tictoc::toc()

# csr (gene-centric) -----------------------------------------------------------

## TODO

csr_matrix <- as(
  Matrix::Matrix(count_mat, sparse = TRUE),
  "RsparseMatrix"
)

length(csc_matrix@p)

csr_matrix@p
csr_matrix@j
csr_matrix@x


rextendr::document()
