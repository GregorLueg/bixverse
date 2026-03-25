# Plot the highly variable genes

Plots the median-absolute deviation of the genes and applied thresholds.
Expects that
[`preprocess_bulk_coexp()`](https://gregorlueg.github.io/bixverse/reference/preprocess_bulk_coexp.md)
was run and will throw an error otherwise.

## Usage

``` r
plot_hvgs(object, bins = 50L)
```

## Arguments

- object:

  The underlying class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

- bins:

  Integer. Number of bins to plot.
