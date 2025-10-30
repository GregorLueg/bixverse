# sc aggregations --------------------------------------------------------------

## testing parameters ----------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
# hvg
hvg_to_keep <- 30L
# pca
no_pcs <- 10L

## synthetic test data ---------------------------------------------------------

single_cell_test_data <- generate_single_cell_test_data()

genes_pass <- which(
  Matrix::colSums(single_cell_test_data$counts != 0) >= min_cells_exp
)

cells_pass <- which(
  (Matrix::rowSums(single_cell_test_data$counts[, genes_pass]) >=
    min_lib_size) &
    (Matrix::rowSums(single_cell_test_data$counts[, genes_pass] != 0) >=
      min_genes_exp)
)

expect_true(
  current = length(genes_pass) > 80 & length(genes_pass) != 100,
  info = "sc processing - sensible amount of genes pass"
)

expect_true(
  current = length(cells_pass) > 800 & length(cells_pass) != 1000,
  info = "sc processing - sensible amount of cells pass"
)

sc_qc_param <- params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

## underlying class ------------------------------------------------------------

sc_object <- single_cell_exp(dir_data = tempdir())

sc_object <- # keep all cells for the sake of this
  sc_object <- load_r_data(
    object = sc_object,
    counts = single_cell_test_data$counts,
    obs = single_cell_test_data$obs,
    var = single_cell_test_data$var,
    sc_qc_param = params_sc_min_quality(
      min_unique_genes = min_genes_exp,
      min_lib_size = min_lib_size,
      min_cells = min_cells_exp
    ),
    streaming = FALSE,
    .verbose = FALSE
  )

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = no_pcs,
  randomised_svd = TRUE,
  .verbose = FALSE
)

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(knn_algorithm = "annoy"),
  .verbose = FALSE
)


object = sc_object
sc_meta_cell_params = params_sc_metacells(target_no_metacells = 50L, k = 5L)
regenerate_knn = FALSE
embd_to_use = "pca"
no_embd_to_use = NULL
target_size = 1e4
seed = 42L
return_aggregated = TRUE
.verbose = TRUE


# if the kNN graph shall be regenerated, get the emedding here...
if (regenerate_knn) {
  embd <- switch(embd_to_use, pca = get_pca_factors(object))
  # early return
  if (is.null(embd)) {
    warning(
      paste(
        "The desired embedding was not found. Please check the parameters.",
        "Returning NULL."
      )
    )

    return(NULL)
  }

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }
  knn_data <- NULL
} else {
  embd <- NULL
  knn_data <- get_knn_mat(object)

  if (is.null(knn_data)) {
    warning(
      paste(
        "No kNN data could be found on the object. Set regenerate_knn to",
        "TRUE or generate the kNN matrix via other means",
        "Returning NULL."
      )
    )
    return(NULL)
  }
}

meta_cell_data <- rs_get_metacells(
  f_path = get_rust_count_cell_f_path(object),
  knn_mat = knn_data,
  embd = embd,
  meta_cell_params = sc_meta_cell_params,
  target_size = target_size,
  seed = seed,
  verbose = .verbose,
  return_aggregated = return_aggregated
)

meta_cell_data$assignments

table(meta_cell_data$assignments$assignments)

meta_cell_data$assignments$metacells


seacell_data <- rs_get_seacells(
  f_path = get_rust_count_cell_f_path(object),
  embd = get_pca_factors(object),
  seacells_params = list(
    n_sea_cells = 50L,
    k = 5L,
    ann_dist = "euclidean",
    convergence_epsilon = 1e-3
  ),
  target_size = target_size,
  seed = seed,
  verbose = .verbose,
  return_aggregated = return_aggregated
)


seacell_data$assignments$assignments

seacell_data$assignments$metacells

table(seacell_data$assignments$assignments)
