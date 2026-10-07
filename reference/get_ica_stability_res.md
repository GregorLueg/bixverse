# Get the ICA component data (stability, convergence, redundancy)

Getter function to extract the ICA component data in terms of stability,
convergence and maximum absolute correlation between the component
loadings (`max_abs_loading_cor`). If not found will return `NULL`.

## Usage

``` r
get_ica_stability_res(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

data.table with the ICA parameter data (if found. Otherwise `NULL`.)

## Examples

``` r
# the stability table behind the ncomp choice
mat <- t(synthetic_signal_matrix()$mat)
obj <- BulkCoExp(mat, data.table::data.table(sample_id = rownames(mat)))
obj <- preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
obj <- ica_processing(obj, .verbose = FALSE)
obj <- ica_evaluate_comp(
  obj,
  ica_type = "logcosh",
  ncomp_params = params_ica_ncomp(custom_seq = seq(2L, 20L, by = 2L)),
  .verbose = FALSE
)
head(get_ica_stability_res(obj))
#> Key: <no_components>
#>    no_components median_stability UQ_stability LQ_stability converged
#>            <int>            <num>        <num>        <num>     <num>
#> 1:             2        0.3341118    0.3371207    0.3311030      1.00
#> 2:             4        0.5896017    0.6422570    0.5311256      1.00
#> 3:             6        0.6287036    0.6477328    0.5436103      0.34
#> 4:             8        0.4903164    0.5704741    0.4185860      0.00
#> 5:            10        0.6038098    0.7605020    0.4277218      0.58
#> 6:            12        0.4654485    0.5977709    0.3903346      0.30
#>    max_abs_loading_cor combined_score
#>                  <num>          <num>
#> 1:          0.01698197     0.32843795
#> 2:          0.62200673     0.22286546
#> 3:          0.75415655     0.05255131
#> 4:          0.70843868     0.00000000
#> 5:          0.84741429     0.05343699
#> 6:          0.84625583     0.02146800
```
