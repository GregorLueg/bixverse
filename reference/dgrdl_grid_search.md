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

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).
  Ideally, you should run
  [`preprocess_bulk_coexp()`](https://gregorlueg.github.io/bixverse/reference/preprocess_bulk_coexp.md)
  before applying this function.

- neighbours_vec:

  Integer vector. The different k nearest neighbours to test.

- dict_size_vec:

  Integer vector. The different dictionary sizes to test for.

- seed_vec:

  Integer vector. The different initial seeds to test for the dictionary
  generation.

- dgrdl_params:

  List. Output of
  [`params_dgrdl()`](https://gregorlueg.github.io/bixverse/reference/params_dgrdl.md):

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

`BulkCoExp` with the grid search results added to the class.

## References

Pan et al., Cell Syst, 2022

## Examples

``` r
# small grid over neighbours and dictionary sizes
syn <- generate_gene_module_data(n_samples = 24L, n_genes = 60L)
obj <- BulkCoExp(syn$data, syn$meta_data)
obj <- preprocess_bulk_coexp(obj, hvg = NULL, .verbose = FALSE)
obj <- dgrdl_grid_search(
  obj,
  neighbours_vec = c(3L, 5L),
  dict_size_vec = c(4L, 6L),
  seed_vec = 123L,
  .verbose = FALSE
)
get_grid_search_res(obj)
#>     seed dict_size k_neighbours reconstruction_errs feature_laplacian_objective
#>    <num>     <num>        <num>               <num>                       <num>
#> 1:   123         4            3            7.809673                  0.08282019
#> 2:   123         4            5            8.081732                  0.06261729
#> 3:   123         6            3            3.477399                  0.09445959
#> 4:   123         6            5            3.794702                  0.07474547
#>    sample_laplacian_objective
#>                         <num>
#> 1:                  0.3560355
#> 2:                  1.7663868
#> 3:                  0.5568773
#> 4:                  2.7064393
```
