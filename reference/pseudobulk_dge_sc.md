# Run the edgeR quasi-likelihood workflow on pseudo-bulked single cells

Sums the raw counts per sample and treats the result as a bulk
experiment. That is the whole method, and it is the one that holds its
nominal false discovery rate when the cells within a sample are not
independent, which they never are.

This is a convenience wrapper: it pseudo-bulks with
[`get_pseudobulked_sc()`](https://gregorlueg.github.io/bixverse/reference/get_pseudobulked_sc.md)
and hands the aggregate to
[`run_edger_ql()`](https://gregorlueg.github.io/bixverse/reference/run_edger_ql.md).
Reach for the two separately if you want to inspect or reuse the
aggregated matrix.

## Usage

``` r
pseudobulk_dge_sc(
  object,
  cell_list,
  design,
  coef = NULL,
  contrast = NULL,
  edger_params = params_edger_ql(),
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- cell_list:

  Named list of character vectors. The cell identifiers per pseudo-bulk
  sample. The names become the sample identifiers, and the rows of
  `design` must follow that order.

- design:

  Numeric matrix. The design matrix of samples x coefficients, including
  the intercept. Rows aligned to `cell_list`.

- coef:

  Optional integer or character. See
  [`run_edger_ql()`](https://gregorlueg.github.io/bixverse/reference/run_edger_ql.md).

- contrast:

  Optional numeric vector or matrix. See
  [`run_edger_ql()`](https://gregorlueg.github.io/bixverse/reference/run_edger_ql.md).

- edger_params:

  A list, see
  [`params_edger_ql()`](https://gregorlueg.github.io/bixverse/reference/params_edger_ql.md).

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A data.table, see
[`run_edger_ql()`](https://gregorlueg.github.io/bixverse/reference/run_edger_ql.md).
The `feature_id` column holds the gene identifiers.

## References

Squair, et al., Nat Commun, 2021; Chen, Lun and Smyth, F1000Research,
2016

## Examples

``` r
# six pseudo-bulk samples out of a demo single cell object
sc <- demo_single_cells(prepped = FALSE)
cells <- get_sc_obs(sc)$cell_id
cell_list <- split(
  cells,
  rep(sprintf("sample_%i", 1:6), length.out = length(cells))
)
design <- stats::model.matrix(~ rep(c("ctrl", "case"), each = 3))
res <- pseudobulk_dge_sc(sc, cell_list, design, .verbose = FALSE)
head(res)
#>    feature_id      log_fc  log_cpm       f_stat   p_value      fdr
#>        <char>       <num>    <num>        <num>     <num>    <num>
#> 1:    gene_01 -0.13045767 14.91690 0.0075305176 0.9331183 0.998572
#> 2:    gene_02 -0.08712794 14.70529 0.0037981562 0.9524644 0.998572
#> 3:    gene_03  0.13829864 14.92486 0.0085170255 0.9288857 0.998572
#> 4:    gene_04  0.09796296 14.50131 0.0048693826 0.9461838 0.998572
#> 5:    gene_05 -0.03489488 13.46428 0.0007964493 0.9782097 0.998572
#> 6:    gene_06 -0.17900623 14.64346 0.0168817897 0.9000273 0.998572

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
