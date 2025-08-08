# sparse data format -----------------------------------------------------------

library(magrittr)

## csc (cell-centric) ----------------------------------------------------------

set.seed(123)

genes <- sprintf("gene_%i", 1:1000)
count_range <- 1:20
no_cells <- 10000L
probs <- 1 / seq_len(length(genes))
probs <- probs / sum(probs)

mirai::daemons(5L)

count_mat <- mirai::mirai_map(
  1:no_cells,
  \(idx, count_range, genes, probs) {
    set.seed(idx)
    random_counts_i <- sample(count_range, 750, replace = TRUE)
    random_indx_i <- sort(sample(1:length(genes), 750, prob = probs))

    counts_i <- rep(0, length(genes))
    counts_i[random_indx_i] <- random_counts_i

    counts_i
  },
  .args = list(count_range, genes, probs)
)[] %>%
  do.call(cbind, .) %>%
  `rownames<-`(genes) %>%
  `colnames<-`(sprintf("cell_%i", 1:no_cells))

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
