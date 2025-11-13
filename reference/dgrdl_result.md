# Run DGRDL with the specified parameters

Runs the DGRDL algorithm from Pan et al., with the specified
hyperparamters. To determine the hyperparameters, you can use
[`dgrdl_grid_search()`](dgrdl_grid_search.md).

## Usage

``` r
dgrdl_result(
  object,
  dgrdl_params = params_dgrdl(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md). Ideally, you should
  run [`preprocess_bulk_coexp()`](preprocess_bulk_coexp.md) before
  applying this function.

- dgrdl_params:

  List. Output of [`params_dgrdl()`](params_dgrdl.md):

  - sparsity - Integer. Sparsity constraint (max non-zero coefficients
    per signal)

  - dict size - Integer. The dictionary size.

  - alpha - Float. Sample context regularisation weight.

  - beta - Float. Feature effect regularisation weight.

  - max_iter - Integer. Maximum number of iterations for the main
    algorithm.

  - k_neighbours - Integer. Number of neighbours for the KNN graph for
    the feature and sample Laplacian.

  - admm_iter - Integer. ADMM iterations for sparse coding.

  - rho - Float. ADMM step size.

- seed:

  Integer. Seed for the initialisation of the dictionary.

- .verbose:

  Boolean. Controls verbosity of the function.

## References

Pan et al., Cell Syst, 2022
