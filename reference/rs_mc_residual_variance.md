# Residual variance and variable features for meta cells

**\[experimental\]** In-memory version of
[`rs_sc_residual_variance()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_residual_variance.md).

## Usage

``` r
rs_mc_residual_variance(sparse_data, residual_fit, n_hvg, verbose)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `cs_type`. Shape is (metacells, genes). Pass raw counts.

- residual_fit:

  List. A fit from
  [`rs_mc_fit_residuals()`](https://gregorlueg.github.io/bixverse/reference/rs_mc_fit_residuals.md).

- n_hvg:

  Integer. Variable features to take from each group.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- genes - The gene indices the variances are indexed by. (0-indexed!)

- variance - Matrix of residual variance, genes by groups.

- hvg - The selected gene indices, ascending. (0-indexed!)
