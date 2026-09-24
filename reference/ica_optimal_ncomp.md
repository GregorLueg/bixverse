# Identify stability inflection point

This function can be used after having run
[`ica_evaluate_comp()`](https://gregorlueg.github.io/bixverse/reference/ica_evaluate_comp.md).
It will calculate the inflection point, based on the first derivative of
a loess function fitted `ncomp ~ combined_score` with the combined score
being a product of the median stability, orthogonality and proportion of
convergence and add these info to the object. Should the loess function
raise a warning (e.g., singularity), the class will be returned as is
and manual determination of optimal ncomp is warranted. Additionally,
you have the option to plot the loess function for additional control
over the span parameter (defaults to `TRUE`).

## Usage

``` r
ica_optimal_ncomp(object, span = 0.2, show_plot = TRUE, .verbose = TRUE)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).
  You need to apply
  [`ica_evaluate_comp()`](https://gregorlueg.github.io/bixverse/reference/ica_evaluate_comp.md)
  before running this function.

- span:

  Float. A value between 0.1 and 1.0 for the loess function defining the
  span. Defaults to `0.2`.

- show_plot:

  Boolean. Shall a plot be shown of the loess function fitted. Defaults
  to `TRUE`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

`BulkCoExp` with optimal ncomp based on the inflection point method.

## Examples

``` r
# inflection point of the combined score
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
obj <- ica_optimal_ncomp(
  obj, span = 0.4, show_plot = FALSE, .verbose = FALSE
)
#> Warning: The loess function could not be fitted with the given parameters. Returning object as is.
obj@params$ica_stability_assessment$optimal_ncomp
#> [1] NA
```
