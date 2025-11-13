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

- dict_size:

  Integer. Dictionary size

- alpha:

  Float. Sample context regularisation weight.

- beta:

  Float. Feature effect regularisation weight.

- max_iter:

  Integer. Maximum number of iterations for the main algorithm.

- k_neighbours:

  Integer. Number of neighbours in the KNN graph.

- admm_iter:

  Integer. ADMM iterations for sparse coding.

- rho:

  Float. ADMM step size.

## Value

List with parameters for usage in subsequent function.
