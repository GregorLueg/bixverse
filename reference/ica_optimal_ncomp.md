# Identify stability inflection point

This function can be used after having run
[`ica_evaluate_comp()`](ica_evaluate_comp.md). It will calculate the
inflection point, based on the first derivative of a loess function
fitted `ncomp ~ combined_score` with the combined score being a product
of the median stability, orthogonality and proportion of convergence and
add these info to the object. Should the loess function raise a warning
(e.g., singularity), the class will be returned as is and manual
determination of optimal ncomp is warrented. Additionally, you have the
option to plot the loess function for additional control over the span
parameter (defaults to `TRUE`).

## Usage

``` r
ica_optimal_ncomp(object, span = 0.2, show_plot = TRUE, .verbose = TRUE)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md). You need to apply
  [`ica_evaluate_comp()`](ica_evaluate_comp.md) before running this
  function.

- span:

  Float. A value between 0.1 and 1.0 for the loess function defining the
  span. Defaults to `0.2`.

- show_plot:

  Boolean. Shall a plot be shown of the loess function fitted. Defaults
  to `TRUE`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

`bulk_coexp` with optimal ncomp based on the inflection point method.
