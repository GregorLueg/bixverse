# Fits a residual model for single cell data

**\[experimental\]** Fits either scTransform (v2) or the analytic
Pearson residual model of Lause, Berens and Kobak over the selected
cells. With `group_of_cell` one model is fitted per group, which is what
a multi-sample experiment wants: each sample keeps its own sequencing
depth and composition.

The fit is returned as a list rather than applied to anything. Hand it
back to
[`rs_sc_residual_variance()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_residual_variance.md),
[`rs_sc_pca_residuals()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_pca_residuals.md)
or
[`rs_sct_corrected_counts()`](https://gregorlueg.github.io/bixverse/reference/rs_sct_corrected_counts.md)
to use it.

## Usage

``` r
rs_sc_fit_residuals(
  f_path_gene,
  f_path_cell,
  method,
  cell_indices,
  group_of_cell,
  covariates,
  params,
  gene_batch_size,
  seed,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- f_path_cell:

  String. Path to the `counts_cells.bin` file. Used for the library
  sizes, and for the per-cell totals of the analytic Pearson fit.

- method:

  String. One of `c("sctransform", "analytic_pearson")`.

- cell_indices:

  Integer vector. The cell indices to use. (0-indexed!)

- group_of_cell:

  Integer vector or `NULL`. Group label per selected cell. (0-indexed,
  dense!) `NULL` fits one model over every cell.

- covariates:

  Named list of numeric vectors, one per covariate, each of length
  `length(cell_indices)`. scTransform only. The order is remembered and
  checked on every subsequent use.

- params:

  Named list. The parameters, see
  [`params_sc_sctransform()`](https://gregorlueg.github.io/bixverse/reference/params_sc_sctransform.md)
  or
  [`params_sc_apr()`](https://gregorlueg.github.io/bixverse/reference/params_sc_apr.md).

- gene_batch_size:

  Integer or `NULL`. Genes held in memory per batch.

- seed:

  Integer. Seed for the step-1 subsample.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- method - The method that was fitted.

- models - The per-group models.

- genes - The gene indices modelled in every group. (0-indexed!)

- group_of_cell - The group label per selected cell. (0-indexed!)

- n_groups - Number of groups.

- cell_indices - The cells that were fitted on. (0-indexed!)

- covariates - The covariates used, scTransform only.

- log10_umi - The per-cell offset, scTransform only.

- cell_totals - The per-cell totals, analytic Pearson only.

## References

Choudhary and Satija, Genome Biology, 2022; Lause, Berens and Kobak,
Genome Biology, 2021.
