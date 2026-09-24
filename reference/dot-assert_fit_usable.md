# Fetch a fitted residual model that is safe to compute with

The getter warns and hands back `NULL` when nothing was fitted, which is
what a presence probe wants. Anything about to compute residuals needs
the harder version, so this errors instead.

Also compares the fitted cell set against the current one. Rust makes
the same check, but by the time it fires the message is about index
vectors; here it can name the function to re-run.

## Usage

``` r
.assert_fit_usable(object, cell_indices)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- cell_indices:

  Integer. The 0-based cells about to be used.

## Value

The `ScResidualFit`.
