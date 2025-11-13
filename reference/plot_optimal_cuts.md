# Plot the k cuts vs median R2

Plots the optimal k vs. median of median R2 graph to identify the
optimal number of cuts.

## Usage

``` r
plot_optimal_cuts(object)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md).

## Value

If optimal cuts results were found, returns the ggplot. Otherwise,
throws a warning and returns NULL.
