# Return the epsilon data

Getter function to extract the
`epsilon param ~ power law goodness of fit` data from the
[`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
class.

## Usage

``` r
get_epsilon_res(object)
```

## Arguments

- object:

  `BulkCoExp` class.

## Value

Returns the epsilon data. (If found. Otherwise `NULL`).

## Examples

``` r
# scale-free fit of the affinity matrix across epsilons
syn <- synthetic_bulk_cor_matrix()
mat <- log1p(t(syn$counts))
meta <- data.table::data.table(sample_id = rownames(mat))
object <- BulkCoExp(raw_data = mat, meta_data = meta)
object <- preprocess_bulk_coexp(object, hvg = 200L, .verbose = FALSE)
object <- cor_module_processing(
  object,
  cor_method = "pearson",
  .verbose = FALSE
)
object <- cor_module_check_epsilon(
  object,
  rbf_func = "gaussian",
  .verbose = FALSE
)
head(get_epsilon_res(object))
#>    epsilon   r2_vals
#>      <num>     <num>
#> 1:    10.0 0.6472736
#> 2:     9.5 0.5733779
#> 3:     9.0 0.6103945
#> 4:     8.5 0.6367963
#> 5:     8.0 0.6010812
#> 6:     7.5 0.4571262
```
