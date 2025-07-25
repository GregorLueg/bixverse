synthetic_data <- generate_gene_module_data(
  n_samples = 500L,
  n_genes = 10000L,
  n_modules = 8L
)


res_simple <- rs_sparse_dict_dgrdl(
  x = synthetic_data,
  dgrdl_params = params_dgrdl(
    sparsity = 10L,
    dict_size = 8L,
    alpha = 1.0,
    beta = 1.0,
    max_iter = 20L,
    k_neighbours = 20L,
    admm_iter = 5L,
    rho = 1.0
  ),
  verbose = TRUE,
  seed = 10101L
)
