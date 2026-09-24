# Generate a sparse dictionary with DGRDL

**\[experimental\]** This is the Rust implementation of dual graph
regularised dictionary learning in the implementation of Pan, et al.,
Cell Systems, 2022. This helper function is designed to run a grid
search over the data.

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

  A list with the parameters for the algorithm. Missing items fall back
  to defaults. Expects the following items.

  - sparsity - Integer. Sparsity constraint (max non-zero coefficients
    per signal).

  - dict_size - Integer. Size of the dictionary. Ignored here,
    `dict_sizes` is used instead.

  - alpha - Float. Sample context regularisation weight. The higher the
    stronger the regularisation.

  - beta - Float. Feature context regularisation weight. The higher the
    stronger the regularisation.

  - max_iter - Integer. Maximum iteration for the algorithm.

  - k_neighbours - Integer. Number of k neighbours for the sample and
    feature Laplacian matrix. Ignored here, `k_neighbours_vec` is used
    instead.

  - admm_iter - Integer. Number of iterations for using alternating
    direction method of multipliers (ADMM).

  - rho - Float. ADMM step size.

- seeds:

  Integer vector. The random seeds to include in the grid search.

- dict_sizes:

  Integer vector. The dictionary sizes to test in the grid search.

- k_neighbours_vec:

  Integer vector. The number of neighbours for the kNN graph generation
  to test in the grid search.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the following elements, one entry per tested combination:

- seed - The tested seeds.

- dict_size - The tested dictionary sizes.

- k_neighbours - The tested numbers of neighbours.

- reconstruction_errs - The reconstruction errors (squared Frobenius
  norm) for these hyperparameters.

- feature_laplacian_objective - The objective values of the feature
  Laplacian term for these hyperparameters.

- sample_laplacian_objective - The objective values of the sample
  Laplacian term for these hyperparameters.
