# Fit a residual model over a set of cells

Shared body of the
[`fit_residuals_sc()`](https://gregorlueg.github.io/bixverse/reference/fit_residuals_sc.md)
methods. Resolves the parameters, the grouping and the covariates, then
hands over to Rust.

## Usage

``` r
.fit_residuals(
  object,
  cell_indices,
  obs,
  method,
  group_column,
  covariate_columns,
  residual_params,
  gene_batch_size,
  seed,
  .verbose
)
```

## Arguments

- object:

  `SingleCells`, `SingleCellsSubset` or `MetaCells` class.

- cell_indices:

  Integer. The 0-based cells to fit on.

- obs:

  data.table. The observation table for those cells, in the same order.

- method:

  String. One of `c("sctransform", "analytic_pearson")`. scTransform
  fits a negative binomial per gene and regularises the parameters;
  analytic Pearson is the closed-form alternative with one shared
  dispersion, which is much cheaper.

- group_column:

  String or `NULL`. Column in the observation table to fit separate
  models over, usually the sample. `NULL` fits one model.

- covariate_columns:

  Character vector or `NULL`. Numeric columns in the observation table
  to add to the design. scTransform only. The library size is never a
  covariate, it enters as a fixed offset.

- residual_params:

  List or `NULL`. Parameters, see
  [`params_sc_sctransform()`](https://gregorlueg.github.io/bixverse/reference/params_sc_sctransform.md)
  or
  [`params_sc_apr()`](https://gregorlueg.github.io/bixverse/reference/params_sc_apr.md).
  `NULL` takes the defaults for `method`.

- gene_batch_size:

  Integer or `NULL`. Genes held in memory per batch.

- seed:

  Integer. Random seed for the step-1 subsample.

- .verbose:

  Boolean or Integer. Controls verbosity.

## Value

The `ScResidualFit`.
