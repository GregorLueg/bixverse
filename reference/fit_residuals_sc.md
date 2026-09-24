# Fit a residual model for single cell data

Fits scTransform (v2) or the analytic Pearson residual model over the
cells currently kept. The fit lands in the cache, and
[`find_hvg_sc()`](https://gregorlueg.github.io/bixverse/reference/find_hvg_sc.md)
with `params_sc_hvg(method = "residual")`,
[`calculate_pca_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_pca_sc.md)
with `residuals = TRUE` and
[`sct_corrected_counts_sc()`](https://gregorlueg.github.io/bixverse/reference/sct_corrected_counts_sc.md)
all read it from there.

Unlike the log-normalised layer, which is computed once at ingestion and
written to disk, nothing is precomputed here: the residual rows are
regenerated on demand by whatever consumes the fit.

With `group_column` one model is fitted per group. That is what a
multi-sample experiment wants, since each sample keeps its own depth and
composition, and it changes how the residual HVG selection behaves, see
[`params_sc_hvg()`](https://gregorlueg.github.io/bixverse/reference/params_sc_hvg.md).

For `MetaCells` the counts are summed UMIs, so the negative binomial
still applies, but at a much greater depth than a single cell. Prefer
`method = "analytic_pearson"` there, and revisit the defaults of
[`params_sc_sctransform()`](https://gregorlueg.github.io/bixverse/reference/params_sc_sctransform.md):
`n_genes` and `n_cells` are sized for raw cells, of which there are
usually far more than meta cells.

## Usage

``` r
fit_residuals_sc(
  object,
  method = c("sctransform", "analytic_pearson"),
  group_column = NULL,
  covariate_columns = NULL,
  residual_params = NULL,
  gene_batch_size = NULL,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`, `SingleCellsSubset` or `MetaCells` class.

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

The object with the fitted model attached.

## References

Choudhary and Satija, Genome Biology, 2022; Lause, Berens and Kobak,
Genome Biology, 2021.

## Examples

``` r
# analytic Pearson residuals over every kept cell
sc <- demo_single_cells(prepped = FALSE)
sc <- fit_residuals_sc(sc, method = "analytic_pearson", .verbose = FALSE)
get_residual_fit(sc)$n_groups
#> [1] 1

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
