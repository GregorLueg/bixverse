synthetic_data <- generate_gene_module_data(
  n_samples = 24L,
  n_genes = 60L,
  n_modules = 4L
)


res_simple <- rs_sparse_dict_dgrdl(
  x = synthetic_data,
  dgrdl_params = params_dgrdl(
    sparsity = 10L,
    dict_size = 8L,
    alpha = 1.0,
    beta = 1.0,
    max_iter = 20L,
    k_neighbours = 5L,
    admm_iter = 5L,
    rho = 1.0
  ),
  verbose = TRUE,
  seed = 10101L
)


neighbours_vector <- as.integer(c(2, 4))
seed_vector <- as.integer(c(123))
dict_size <- as.integer(c(4, 6, 8))

tictoc::tic()
grid_search_res <- rs_sparse_dict_dgrdl_grid_search(
  x = synthetic_data,
  dgrdl_params = params_dgrdl(
    sparsity = 10L,
    dict_size = 8L,
    alpha = 1.0,
    beta = 1.0,
    max_iter = 10L,
    k_neighbours = 5L,
    admm_iter = 5L,
    rho = 1.0
  ),
  seeds = seed_vector,
  dict_sizes = dict_size,
  k_neighbours_vec = neighbours_vector,
  verbose = FALSE
)
tictoc::toc()

grid_search_res_dt <- as.data.table(grid_search_res)

ggplot(
  data = grid_search_res_dt,
  mapping = aes(
    x = k_neighbours,
    y = feature_laplacian_objective,
    size = dict_size
  )
) +
  geom_point()
