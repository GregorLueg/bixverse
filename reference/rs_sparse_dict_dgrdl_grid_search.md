# Generate a sparse dictionary with DGRDL

This is the Rust implementation of dual graph regularised dictionary
learning in the implementation of Pan, et al., Cell Systems, 2022. This
helper function is designed to run a grid search over the data.

## Usage

``` r
rs_sparse_dict_dgrdl_grid_search(
  x,
  dgrdl_params,
  seeds,
  dict_sizes,
  k_neighbours_vec,
  verbose
)
```

## Arguments

- x:

  Numerical matrix. Rows = samples, columns = features.

- dgrdl_params:

  A list with the parameters for the algorithm. Expects the following
  items.

  - sparsity - Sparsity constraint (max non-zero coefficients per
    signal).

  - dict_size - Size of the dictionary. This parameter will be ignored
    for this function and `dict_sizes` will be used.

  - alpha - Float. Sample context regularisation weight. The higher the
    stronger the regularisation.

  - beta - Float. Feature context regularisation weight. The higher the
    stronger the regularisation.

  - max_iter - Integer. Maximum iteration for the algorithm.

  - k_neighbours - Integer. Number of k neighbours for the sample and
    feature Laplacian matrix for the regularisation. This parameter will
    be ignored and `k_neighbours_vec` will be used.

  - admm_iter Integer. Number of iterations for using alternating
    direction method of multipliers (ADMM).

  - rho Float. ADMM step size.

- seeds:

  Integer vectors. The random seeds to include in the grid search.

- dict_sizes:

  Integer vector. The dictionary sizes to test in the grid search.

- k_neighbours_vec:

  Integer vector. The number of neighbours for the KNN graph generation
  to test in the grid search.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the following elements:

- seed - The tested seeds.

- dict_size - The tested dictionary sizes.

- reconstruction_errs - The reconstruction errors for these hyper
  parameters.

- feature_laplacian_objective - The objective values of the feature
  Laplacian term for these hyperparameters.

- sample_laplacian_objective - The objective values of the sample
  Laplacian term for these hyperparameters.
