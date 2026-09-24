# Wrapper function to generate DGRDL parameters

Wrapper function to generate DGRDL parameters

## Usage

``` r
params_dgrdl(
  sparsity = 5L,
  dict_size = 5L,
  alpha = 1,
  beta = 1,
  max_iter = 20L,
  k_neighbours = 5L,
  admm_iter = 5L,
  rho = 1
)
```

## Arguments

- sparsity:

  Integer. Sparsity constraint (max non-zero coefficients per signal)
  Defaults to `5L`.

- dict_size:

  Integer. Dictionary size Defaults to `5L`.

- alpha:

  Numeric. Sample context regularisation weight. Defaults to `1.0`.

- beta:

  Numeric. Feature effect regularisation weight. Defaults to `1.0`.

- max_iter:

  Integer. Maximum number of iterations for the main algorithm. Defaults
  to `20L`.

- k_neighbours:

  Integer. Number of neighbours in the KNN graph. Defaults to `5L`.

- admm_iter:

  Integer. ADMM iterations for sparse coding. Defaults to `5L`.

- rho:

  Numeric. ADMM step size. Defaults to `1.0`.

## Value

A named list with the following elements:

- sparsity - Integer. Sparsity constraint (max non-zero coefficients per
  signal) Defaults to `5L`.

- dict_size - Integer. Dictionary size Defaults to `5L`.

- alpha - Numeric. Sample context regularisation weight. Defaults to
  `1.0`.

- beta - Numeric. Feature effect regularisation weight. Defaults to
  `1.0`.

- max_iter - Integer. Maximum number of iterations for the main
  algorithm. Defaults to `20L`.

- k_neighbours - Integer. Number of neighbours in the KNN graph.
  Defaults to `5L`.

- admm_iter - Integer. ADMM iterations for sparse coding. Defaults to
  `5L`.

- rho - Numeric. ADMM step size. Defaults to `1.0`.
