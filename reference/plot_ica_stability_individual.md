# Plot the stability of the ICA components

Helper function to plot the individual stability profiles over the
tested ncomps.

## Usage

``` r
plot_ica_stability_individual(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).
  You need to apply
  [`ica_evaluate_comp()`](https://gregorlueg.github.io/bixverse/reference/ica_evaluate_comp.md)
  before running this function.

## Value

A ggplot with the per-component stability profiles over the tested
number of components.

## Examples

``` r
# per-component stability profiles across the tested ncomps
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
plot_ica_stability_individual(obj)
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_line()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_point()`).
```
