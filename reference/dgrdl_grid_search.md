# Grid search over DGRDL parameters

This function allows you to quickly iterate over different initial
seeds, number of neighbours for the KNN graph and dictionary sizes to
identify optimal hyperparameters for your DGRDL run.

## Usage

``` r
dgrdl_grid_search(
  object,
  neighbours_vec,
  dict_size_vec,
  seed_vec,
  dgrdl_params = params_dgrdl(),
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md). Ideally, you should
  run [`preprocess_bulk_coexp()`](preprocess_bulk_coexp.md) before
  applying this function.

- neighbours_vec:

  Integer vector. The different k nearest neighbours to test.

- dict_size_vec:

  Integer vector. The different dictionary sizes to test for.

- seed_vec:

  Integer vector. The different initial seeds to test for the dictionary
  generation.

- dgrdl_params:

  List. Output of [`params_dgrdl()`](params_dgrdl.md):

  - sparsity - Integer. Sparsity constraint (max non-zero coefficients
    per signal)

  - dict size - Integer. Will be ignored by this function and the
    `dict_size_vec` vector will be used.

  - alpha - Float. Sample context regularisation weight.

  - beta - Float. Feature effect regularisation weight.

  - max_iter - Integer. Maximum number of iterations for the main
    algorithm.

  - k_neighbours - Integer. Will be ignored by this function and the
    `neighbours_vec` will be used.

  - admm_iter - Integer. ADMM iterations for sparse coding.

  - rho - Float. ADMM step size.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

`bulk_coexp` with the grid search results added to the class.

## References

Pan et al., Cell Syst, 2022
