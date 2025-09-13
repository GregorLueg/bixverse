# sparse data format -----------------------------------------------------------

library(magrittr)

rextendr::document()

## csr (cell-centric) ----------------------------------------------------------

# based on the h5ad format with cells -> rows; genes -> columns

# let's try with 1m cells...
seed <- 123L
no_genes <- 20000L
no_cells <- 1000000L

rs_sparse_data <- rs_synthetic_sc_data_csr(
  n_genes = no_genes,
  n_cells = no_cells,
  min_genes = 250,
  max_genes = 1000,
  max_exp = 50,
  seed = seed
)

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

# rextendr::document()

tictoc::tic()
single_cell_counts$r_csr_mat_to_file(
  no_cells = no_cells,
  no_genes = no_genes,
  data = as.integer(rs_sparse_data$data),
  row_ptr = as.integer(rs_sparse_data$indptr),
  col_idx = as.integer(rs_sparse_data$indices),
  target_size = 1e5
)
tictoc::toc()

tictoc::tic()
single_cell_counts$generate_gene_based_data()
tictoc::toc()

indices <- sort(sample(1:no_cells, 100000))

tictoc::tic()
return_data <- single_cell_counts$get_cells_by_indices(
  indices = order(indices),
  assay = "norm"
)
tictoc::toc()

file.size(f_path_cells) / 1024^2
file.size(f_path_genes) / 1024^2

# csc (gene-centric) -----------------------------------------------------------

gene_indices <- sample(1:no_genes, 1000)

tictoc::tic()
return_gene_data <- single_cell_counts$get_genes_by_indices(
  indices = gene_indices,
  assay = "raw"
)
tictoc::toc()

return_gene_data$row_ptr

# h5 files ---------------------------------------------------------------------

library(duckdb)

h5_path <- "~/Downloads/ERX11148735.h5ad"

get_h5ad_dimensions(h5_path)

list.files(tempdir())

db_path <- file.path(tempdir(), "exp.db")
single_cell_db <- dbConnect(duckdb::duckdb(), dbdir = db_path)


obs <- h5_content[
  group == "/obs" & otype == "H5I_DATASET"
][1, as.numeric(dim)]

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


counts


file.exists("~/Downloads/ERX11148735.h5ad")


rextendr::document()
expanded_path <- path.expand("~/Downloads/ERX11148735.h5ad")

get_h5ad_dimensions(expanded_path)

h5_content <- rhdf5::h5ls(
  expanded_path
) %>%
  data.table::setDT()

counts <- h5_content[
  group == "/X" & name == "indptr"
]

tictoc::tic()
data <- rhdf5::h5read(
  file = expanded_path,
  name = "/X/data"
)

indices <- rhdf5::h5read(
  file = expanded_path,
  name = "/X/indices"
)

indptr <- rhdf5::h5read(
  file = expanded_path,
  name = "/X/indptr"
)
tictoc::toc()

tictoc::tic()
h5ad_data <- rs_h5ad_data(expanded_path, "CSR")
tictoc::toc()

length(h5ad_data$indptr)

get_h5ad_dimensions(expanded_path)

rextendr::document()

devtools::load_all()

db_connection <- single_cell_duckdb_con$new(db_dir = tempdir())

db_connection$populate_obs_from_h5(h5_path = h5_path)$populate_vars_from_h5(
  h5_path = h5_path
)

db_connection$get_obs_table()

db_connection$get_vars_table()
