# Wrapper function for parameters for the edgeR quasi-likelihood workflow

Parameters for the edgeR quasi-likelihood chain, implemented in Rust via
the `edge-rs` crate and gated against edgeR 4.8.2. Defaults are edgeR's
own.

The `legacy` switch picks between two genuinely different pipelines. The
current route estimates its own dispersion from the most abundant genes
and skips `estimateDisp()`, which is where most of the runtime went and
is edgeR 4's own recommendation. The legacy route shrinks the raw
residual deviance, needs a dispersion handed to it, and is the only one
where the Poisson bound bites.

## Usage

``` r
params_edger_ql(
  norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
  filter = TRUE,
  min_mean = 0,
  robust = FALSE,
  legacy = FALSE
)
```

## Arguments

- norm_method:

  String. Library size normalisation. One of
  `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`. Defaults to
  `"TMM"`. `"none"` leaves every factor at one, which is what Milo's
  `logMS` amounts to.

- filter:

  Boolean. Run `filterByExpr()` before fitting. Defaults to `TRUE`. Turn
  this off for anything that is not gene expression, e.g. Milo
  neighbourhood counts, where the heuristic means nothing.

- min_mean:

  Numeric. Drop features whose mean count across samples is below this.
  Applied on top of `filter`. Defaults to `0` (off).

- robust:

  Boolean. Robust empirical Bayes squeezing, giving outlier features
  their own smaller prior degrees of freedom. Defaults to `FALSE`.

- legacy:

  Boolean. Take edgeR's pre-4.0 quasi-likelihood pipeline. Defaults to
  `FALSE`.

## Value

A list with the edgeR quasi-likelihood parameters.

## References

Chen, Lun and Smyth, F1000Research, 2016
