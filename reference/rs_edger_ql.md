# Run the edgeR quasi-likelihood chain on a count matrix

**\[experimental\]** Runs `filterByExpr` -\> `calcNormFactors` -\>
`glmQLFit` -\> `glmQLFTest`, implemented in Rust via the `edge-rs` crate
and gated against edgeR 4.8.2. The tested axis does not have to be
genes: Milo's neighbourhood counts are tested with the same call, with
`filter = FALSE`.

## Usage

``` r
rs_edger_ql(counts, design, edger_params)
```

## Arguments

- counts:

  Numeric matrix. Raw counts of features x samples. Must not be
  normalised or log-transformed.

- design:

  Numeric matrix. The design matrix of samples x coefficients, including
  the intercept. Needs at least two columns, since the null model has to
  retain one.

- edger_params:

  Named list. The edgeR parameters, see
  [`params_edger_ql()`](https://gregorlueg.github.io/bixverse/reference/params_edger_ql.md),
  plus either `coef` (0-indexed(!) design columns to drop from the null
  model) or `contrast` (column-major weights with `n_contrasts`
  columns).

## Value

A list with the following elements

- features_to_keep - Boolean. Which features survived the filters. Spans
  the full feature axis of `counts`.

- log_fc - Log2 fold changes of the tested coefficient or contrast.

- log_cpm - Average log2 counts per million.

- f_stat - The quasi-likelihood F statistic.

- p_values - Raw p-values.

- fdr - Benjamini-Hochberg adjusted p-values.

## References

Chen, Lun and Smyth, F1000Research, 2016
