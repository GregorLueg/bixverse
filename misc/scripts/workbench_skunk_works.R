# metacells and seacell --------------------------------------------------------

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

sc_object[[]]

target_size = 1e4
seed = 42L
.verbose = TRUE
return_aggregated = TRUE

rextendr::document()

tictoc::tic()
seacell_data <- rs_get_seacells(
  f_path = get_rust_count_cell_f_path(sc_object),
  embd = get_pca_factors(sc_object),
  seacells_params = list(
    n_sea_cells = 90L,
    k = 15L,
    convergence_epsilon = 1e-5,
    ann_dist = "euclidean"
  ),
  target_size = target_size,
  seed = seed,
  verbose = .verbose,
  return_aggregated = return_aggregated
)
tictoc::toc()

# results_pre_refactor_rss <- seacell_data$rss

# results_pre_refactor_table <- sort(
#   table(seacell_data$assignments$assignments),
#   decreasing = TRUE
# )

seacell_data$assignments$assignments[1:20]
