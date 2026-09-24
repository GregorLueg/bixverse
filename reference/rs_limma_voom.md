# Run the limma linear model chain on a count matrix

**\[experimental\]** Runs optional `filterByExpr` -\> `calcNormFactors`
-\> `voomLmFit` (or limma-trend) -\> `contrasts.fit` -\> `eBayes` -\>
`topTable` for one coefficient or contrast, implemented in Rust via the
`edge-rs` crate and gated against limma 3.66.0.

## Usage

``` r
rs_limma_voom(counts, design, lib_size, limma_params)
```

## Arguments

- counts:

  Integer or double matrix. Raw counts of genes x samples. Must not be
  normalised or log-transformed.

- design:

  Numeric matrix. The design matrix of samples x coefficients, including
  the intercept. Must be full rank.

- lib_size:

  Numeric vector or NULL. Library size per sample. NULL uses the column
  sums of `counts`. Pass the column sums from before gene filtering to
  match edgeR, which keeps those on a subset `DGEList`.

- limma_params:

  Named list. The limma parameters, see
  [`params_limma_voom()`](https://gregorlueg.github.io/bixverse/reference/params_limma_voom.md),
  plus either `coef` (a single 0-indexed(!) design column) or `contrast`
  (column-major weights with `n_contrasts` columns).

## Value

A list with the following elements, all but `features_to_keep` with one
entry per kept gene, in input order

- features_to_keep - Boolean. Which genes survived the filters. Spans
  the full gene axis of `counts`.

- log_fc - Log2 fold changes of the tested coefficient or contrast.

- ci_lower - Lower end of the 95% confidence interval on `log_fc`.

- ci_upper - Upper end of the 95% confidence interval on `log_fc`.

- ave_expr - Average log2 counts per million.

- t_stat - Moderated t statistic.

- p_values - Raw p-values.

- fdr - Benjamini-Hochberg adjusted p-values.

- b_stat - Log-odds of differential expression.

## References

Law, et al., Genome Biol, 2014; Smyth, Stat Appl Genet Mol Biol, 2004
