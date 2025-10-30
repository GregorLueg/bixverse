# test seacell -----------------------------------------------------------------

data_dir <- path.expand("~/Desktop/seacell_test/")

sc_object <- single_cell_exp(dir_data = data_dir)

h5_path <- path.expand(
  "~/repos/other/SEACells/notebooks/data/cd34_multiome_rna.h5ad"
)

sc_object <- load_h5ad(
  sc_object,
  h5_path = path.expand(
    "~/repos/other/SEACells/notebooks/data/cd34_multiome_rna.h5ad"
  ),
  streaming = FALSE
)

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = 1500L,
  .verbose = FALSE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 50L,
  randomised_svd = TRUE,
  .verbose = FALSE
)

rextendr::document()

sea_cells <- get_seacells_sc(
  sc_object,
  seacell_params = params_sc_seacells(
    n_sea_cells = 90L,
    k = 15L,
    convergence_epsilon = 1e-5,
    max_fw_iters = 50L
  )
)
