# Run the edgeR quasi-likelihood workflow

Runs `filterByExpr()` -\> `calcNormFactors()` -\> `glmQLFit()` -\>
`glmQLFTest()` in Rust via the `edge-rs` crate, gated against edgeR
4.8.2.

The tested axis does not have to be genes. Anything with a counts matrix
of features x samples goes through here, which is why the Milo
neighbourhood test calls the same function with `filter = FALSE`.

By default this skips `estimateDisp()` and lets the fit estimate its own
dispersion from the most abundant features. That is edgeR 4's own
recommendation and where most of the runtime went. Set `legacy = TRUE`
in
[`params_edger_ql()`](https://gregorlueg.github.io/bixverse/reference/params_edger_ql.md)
for the pre-4.0 pipeline, which does run `estimateDisp()`.

## Usage

``` r
run_edger_ql(
  counts,
  design,
  coef = NULL,
  contrast = NULL,
  edger_params = params_edger_ql()
)
```

## Arguments

- counts:

  Numeric matrix. Raw counts of features x samples, with rownames. Must
  not be normalised or log-transformed.

- design:

  Numeric matrix. The design matrix of samples x coefficients, including
  the intercept. Usually the output of
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).
  Needs at least two columns, since the null model has to retain one.

- coef:

  Optional integer or character. Which coefficient(s) of `design` to
  drop from the null model, given as 1-based column positions or column
  names. Defaults to the last column, as edgeR does.

- contrast:

  Optional numeric vector or matrix. Weights over the coefficients, one
  entry (or row) per column of `design`. Mutually exclusive with `coef`.

- edger_params:

  A list, see
  [`params_edger_ql()`](https://gregorlueg.github.io/bixverse/reference/params_edger_ql.md).
  The list has the following parameters:

  - norm_method - String. Library size normalisation. One of
    `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`.

  - filter - Boolean. Run `filterByExpr()` before fitting.

  - min_mean - Numeric. Minimum mean count across samples.

  - robust - Boolean. Robust empirical Bayes squeezing.

  - legacy - Boolean. Take edgeR's pre-4.0 pipeline.

## Value

A data.table with one row per feature that survived the filters:

- feature_id - The rowname of `counts`.

- log_fc - Log2 fold change of the tested coefficient or contrast.

- log_cpm - Average log2 counts per million.

- f_stat - The quasi-likelihood F statistic.

- p_value - Raw p-value.

- fdr - Benjamini-Hochberg adjusted p-value.

## References

Chen, Lun and Smyth, F1000Research, 2016
