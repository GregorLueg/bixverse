# Wrapper function for parameters for the limma-voom workflow

Parameters for the limma linear model chain, implemented in Rust via the
`edge-rs` crate and gated against limma 3.66.0. Defaults are limma's
own, except for `filter`: inside
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md)
the genes were already filtered by
[`qc_bulk_dge()`](https://gregorlueg.github.io/bixverse/reference/qc_bulk_dge.md).
`route = "voom"` is `voomLmFit()`: precision weights from the
mean-variance trend, then weighted least squares. `route = "trend"` is
limma-trend: log-CPM straight into `lmFit()`, with the trend absorbed by
`eBayes(trend = TRUE)`. The empirical Bayes trend follows the route.

## Usage

``` r
params_limma_voom(
  route = c("voom", "trend"),
  norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
  filter = FALSE,
  min_mean = 0,
  robust = FALSE,
  prior_count = NULL,
  adaptive_span = TRUE,
  span = 0.5,
  proportion = 0.01
)
```

## Arguments

- route:

  String. Whether to run the voom or the limma-trend route. One of
  `c("voom", "trend")`. Defaults to `"voom"`.

- norm_method:

  String. Library size normalisation. One of
  `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`. Defaults to
  `"TMM"`.

- filter:

  Boolean. Run `filterByExpr()` before fitting. Defaults to `FALSE`.

- min_mean:

  Numeric. Drop genes whose mean count across samples is below this.
  Applied on top of `filter`. Defaults to `0.0`.

- robust:

  Boolean. Robust empirical Bayes, `eBayes(robust = TRUE)`. Defaults to
  `FALSE`.

- prior_count:

  Numeric or `NULL`. Count added before the log. `NULL` takes the
  route's own default, `0.5` for voom and `2` for trend. Defaults to
  `NULL`.

- adaptive_span:

  Boolean. Derive the lowess span from the number of genes, as limma
  does since 3.56. Only used by voom. Defaults to `TRUE`.

- span:

  Numeric. Lowess span for the voom trend, only read if
  `adaptive_span = FALSE`. Defaults to `0.5`.

- proportion:

  Numeric. Assumed proportion of differentially expressed genes, only
  used for the B-statistic. Defaults to `0.01`.

## Value

A named list with the following elements:

- route - String. Whether to run the voom or the limma-trend route. One
  of `c("voom", "trend")`. Defaults to `"voom"`.

- norm_method - String. Library size normalisation. One of
  `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`. Defaults to
  `"TMM"`.

- filter - Boolean. Run `filterByExpr()` before fitting. Defaults to
  `FALSE`.

- min_mean - Numeric. Drop genes whose mean count across samples is
  below this. Applied on top of `filter`. Defaults to `0.0`.

- robust - Boolean. Robust empirical Bayes, `eBayes(robust = TRUE)`.
  Defaults to `FALSE`.

- prior_count - Numeric or `NULL`. Count added before the log. `NULL`
  takes the route's own default, `0.5` for voom and `2` for trend.
  Defaults to `NULL`.

- adaptive_span - Boolean. Derive the lowess span from the number of
  genes, as limma does since 3.56. Only used by voom. Defaults to
  `TRUE`.

- span - Numeric. Lowess span for the voom trend, only read if
  `adaptive_span = FALSE`. Defaults to `0.5`.

- proportion - Numeric. Assumed proportion of differentially expressed
  genes, only used for the B-statistic. Defaults to `0.01`.

## References

Law, et al., Genome Biol, 2014; Smyth, Stat Appl Genet Mol Biol, 2004
