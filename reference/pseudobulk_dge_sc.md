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
