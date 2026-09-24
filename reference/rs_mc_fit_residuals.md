# Fits a residual model for meta cells

**\[experimental\]** In-memory version of
[`rs_sc_fit_residuals()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_fit_residuals.md).
Meta cell counts are summed UMIs, so the negative binomial the residual
models describe is still defined; it just sits at a much greater depth
than a single cell.

## Usage

``` r
rs_mc_fit_residuals(
  sparse_data,
  method,
  group_of_cell,
  covariates,
  params,
  seed,
  verbose
)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `cs_type`. Shape is (metacells, genes). Pass raw counts.

- method:

  String. One of `c("sctransform", "analytic_pearson")`.

- group_of_cell:

  Integer vector or `NULL`. Group label per meta cell. (0-indexed,
  dense!) `NULL` fits one model.

- covariates:

  Named list of numeric vectors, one per covariate, each of length
  `nrow`. scTransform only.

- params:

  Named list. See
  [`params_sc_sctransform()`](https://gregorlueg.github.io/bixverse/reference/params_sc_sctransform.md)
  or
  [`params_sc_apr()`](https://gregorlueg.github.io/bixverse/reference/params_sc_apr.md).

- seed:

  Integer. Seed for the step-1 subsample.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list as described in
[`rs_sc_fit_residuals()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_fit_residuals.md).
