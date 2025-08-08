# sparse data format -----------------------------------------------------------

library(magrittr)

## csc (cell-centric) ----------------------------------------------------------

set.seed(123)

genes <- letters
count_range <- 1:20
no_cells <- 25L
probs <- 1 / seq_len(length(genes))
probs <- probs / sum(probs)

count_mat <- purrr::map(1:no_cells, \(idx) {
  set.seed(idx)
  random_counts_i <- sample(count_range, 5, replace = TRUE)
  random_indx_i <- sort(sample(1:length(genes), 5, prob = probs))

  counts_i <- rep(0, length(genes))
  counts_i[random_indx_i] <- random_counts_i

  counts_i
}) %>%
  do.call(cbind, .) %>%
  `rownames<-`(letters) %>%
  `colnames<-`(sprintf("cell_%i", 1:no_cells))

storage.mode(count_mat) <- "integer"

csc_matrix <- as(
  Matrix::Matrix(count_mat, sparse = TRUE),
  "CsparseMatrix"
)

storage.mode(csc_matrix)

rextendr::document()

csc_matrix@i
csc_matrix@p
csc_matrix@x
csc_matrix@Dim

length(csc_matrix@p)

dir <- tempdir()
f_path <- file.path(dir, "test.bin")

rs_csc_to_binary_f(
  f_path = f_path,
  no_cells = csc_matrix@Dim[2],
  no_genes = csc_matrix@Dim[1],
  data = as.integer(csc_matrix@x),
  col_ptr = csc_matrix@p,
  row_idx = csc_matrix@i
)

read_file <- rs_binary_f_to_csc(f_path)

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
