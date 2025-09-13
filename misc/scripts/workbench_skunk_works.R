# sparse data format -----------------------------------------------------------

library(magrittr)

rextendr::document()

## csr (cell-centric) ----------------------------------------------------------

# based on the h5ad format with cells -> rows; genes -> columns

# let's try with 1m cells...
seed <- 123L
no_genes <- 20000L
no_cells <- 100000L

tictoc::tic()
rs_sparse_data <- rs_synthetic_sc_data_csr(
  n_genes = no_genes,
  n_cells = no_cells,
  min_genes = 250,
  max_genes = 1000,
  max_exp = 50,
  seed = seed
)
tictoc::toc()

gc()

# csr_matrix <- Matrix::sparseMatrix(
#   j = rs_sparse_data$col_indices + 1,
#   p = rs_sparse_data$row_ptrs,
#   x = rs_sparse_data$data,
#   dims = c(no_cells, no_genes),
#   dimnames = list(
#     sprintf("cell_%i", 1:no_cells),
#     sprintf("gene_%i", 1:no_genes)
#   )
# )

# rm(rs_sparse_data)

# csr_matrix <- as(csr_matrix, "RsparseMatrix")

dir <- tempdir()
f_path_cells <- file.path(dir, "cells.bin")
f_path_genes <- file.path(dir, "genes.bin")

single_cell_counts <- SingeCellCountData$new(
  f_path_cells = f_path_cells,
  f_path_genes = f_path_genes
)

list.files(dir)

tictoc::tic()
single_cell_counts$r_csr_mat_to_file(
  no_cells = no_cells,
  no_genes = no_genes,
  data = as.integer(rs_sparse_data$data),
  row_ptr = as.integer(rs_sparse_data$row_ptrs),
  col_idx = as.integer(rs_sparse_data$col_indices),
  target_size = 1e5
)
tictoc::toc()

tictoc::tic()
single_cell_counts$generate_gene_based_data()
tictoc::toc()

indices <- sample(1:no_cells, 10000)

tictoc::tic()
return_data <- single_cell_counts$get_cells_by_indices(
  indices = order(indices),
  assay = "norm"
)
tictoc::toc()

file.size(f_path_cells) / 1024^2
file.size(f_path_genes) / 1024^2

# csc (gene-centric) -----------------------------------------------------------

gene_indices <- sample(1:no_genes, 100)

return_gene_data <- single_cell_counts$get_genes_by_indices(
  indices = gene_indices,
  assay = "raw"
)

return_gene_data$row_ptr

# h5 files ---------------------------------------------------------------------

library(duckdb)

h5_path <- "~/Downloads/anndata.h5ad"

list.files(tempdir())

db_path <- file.path(tempdir(), "exp.db")
single_cell_db <- dbConnect(duckdb::duckdb(), dbdir = db_path)

h5_content <- rhdf5::h5ls(
  h5_path
) %>%
  data.table::setDT()

obs <- h5_content[
  group == "/obs" & otype == "H5I_DATASET",
  setNames((paste(group, name, sep = "/")), name)
]

names(obs) <- gsub("-", "_", names(obs))

for (i in seq_along(obs)) {
  col_name <- names(obs)[i]
  col_path <- obs[[i]]
  col_data <- data.table(x = rhdf5::h5read(h5_path, col_path)) %>%
    `names<-`(col_name) %>%
    .[, cell_idx := .I] %>%
    .[, c("cell_idx", col_name), with = FALSE]

  if (i == 1) {
    DBI::dbWriteTable(
      single_cell_db,
      "obs",
      col_data,
      overwrite = TRUE
    )
  } else {
    DBI::dbWriteTable(single_cell_db, "temp_col", col_data, overwrite = TRUE)

    DBI::dbExecute(
      single_cell_db,
      sprintf(
        'CREATE TABLE obs_new AS
        SELECT obs.*, temp_col.%s
        FROM obs
        JOIN temp_col ON obs.cell_idx = temp_col.cell_idx;
        DROP table obs;
        ALTER TABLE obs_new RENAME TO obs;
        DROP TABLE temp_col',
        col_name
      )
    )
  }
}

obs <- data.table::setDT(DBI::dbGetQuery(single_cell_db, 'SELECT * FROM obs'))

DBI::dbDisconnect(single_cell_db)

for (i in seq_along(obs)) {
  col_name <- names(obs)[i]
  col_path <- obs[[i]]
  col_data <- data.table(x = rhdf5::h5read(h5_path, col_path)) %>%
    `names<-`(col_name)

  if (i == 1) {
    DBI::dbWriteTable(con, "obs", col_data)
  } else {}

  if (first_column) {
    DBI::dbExecute(
      con,
      sprintf(
        'COPY (SELECT * FROM col_data) TO "%s" (FORMAT PARQUET)',
        parquet_path
      )
    )
    first_column <- FALSE
  } else {}
}


counts <- h5_content[
  group == "/X",
  setNames((paste(group, name, sep = "/")), name)
]


counts

data <- rhdf5::h5read(
  file = h5_path,
  name = counts["data"]
)

indices <- rhdf5::h5read(
  file = h5_path,
  name = counts["indices"]
)

length(indptr)

max(indices)

indices[3618]

indptr <- rhdf5::h5read(
  file = h5_path,
  name = counts["indptr"]
)
