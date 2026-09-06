# Run DGRDL with the specified parameters

Runs the DGRDL algorithm from Pan et al., with the specified
hyperparamters. To determine the hyperparameters, you can use
[`dgrdl_grid_search()`](https://gregorlueg.github.io/bixverse/reference/dgrdl_grid_search.md).

## Usage

``` r
dgrdl_result(
  object,
  dgrdl_params = params_dgrdl(),
  membership_params = params_module_membership(),
  seed = 42L,
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

- dgrdl_params:

  List. Output of
  [`params_dgrdl()`](https://gregorlueg.github.io/bixverse/reference/params_dgrdl.md):

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

- membership_params:

  List. Controls how the atom loadings are turned into module
  membership, see
  [`params_module_membership()`](https://gregorlueg.github.io/bixverse/reference/params_module_membership.md).
  Membership is not exclusive: a gene active in several atoms appears in
  several modules, and a gene in no tail appears in none.

- seed:

  Integer. Seed for the initialisation of the dictionary.

- .verbose:

  Boolean. Controls verbosity of the function.

## References

Pan et al., Cell Syst, 2022

## Examples

``` r
# fit DGRDL with a six-atom dictionary
syn <- generate_gene_module_data(n_samples = 24L, n_genes = 60L)
obj <- BulkCoExp(syn$data, syn$meta_data)
obj <- preprocess_bulk_coexp(obj, hvg = NULL, .verbose = FALSE)
obj <- dgrdl_result(
  obj,
  dgrdl_params = params_dgrdl(dict_size = 6L, k_neighbours = 3L),
  .verbose = FALSE
)
head(get_modules(get_results(obj)))
#>          gene module_id  loading   sign        z
#>        <char>    <char>    <num> <char>    <num>
#> 1: feature_27    dict_2 1.889059    pos 3.077097
#> 2: feature_24    dict_2 1.885384    pos 3.070073
#> 3: feature_22    dict_2 1.882083    pos 3.063763
#> 4: feature_28    dict_2 1.881407    pos 3.062472
#> 5: feature_25    dict_2 1.878301    pos 3.056535
#> 6: feature_23    dict_2 1.877537    pos 3.055074
```
