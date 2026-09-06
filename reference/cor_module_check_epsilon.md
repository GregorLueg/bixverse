# Iterate through different epsilon parameters

This functions iterates through a set of provided epsilons and checks
for each one to what extend the underlying affinity matrix will follow a
power law distribution.

## Usage

``` r
cor_module_check_epsilon(
  object,
  rbf_func = c("bump", "gaussian", "inverse_quadratic"),
  epsilons = c(0.25, seq(from = 0.5, to = 10, by = 0.5)),
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).
  You need to run
  [`cor_module_processing()`](https://gregorlueg.github.io/bixverse/reference/cor_module_processing.md)
  before running this function.

- rbf_func:

  The type of RBF function you want to apply. A choice of
  `c('bump', 'gaussian', 'inverse_quadratic')`.

- epsilons:

  Vector of floats. The different epsilon parameters you would like to
  run.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added data to the properties for subsequent usage.

## Examples

``` r
# scan epsilons of a Gaussian RBF for scale-free topology
mat <- t(synthetic_signal_matrix()$mat)
obj <- BulkCoExp(mat, data.table::data.table(sample_id = rownames(mat)))
obj <- preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
obj <- cor_module_processing(obj, cor_method = "spearman", .verbose = FALSE)
obj <- cor_module_check_epsilon(
  obj, rbf_func = "gaussian", .verbose = FALSE
)
head(get_epsilon_res(obj))
#>    epsilon   r2_vals
#>      <num>     <num>
#> 1:    10.0 0.4062315
#> 2:     9.5 0.4551064
#> 3:     9.0 0.5708543
#> 4:     8.5 0.6176994
#> 5:     8.0 0.6518912
#> 6:     7.5 0.6920401
```
