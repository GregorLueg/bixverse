# Wrapper function to generate module membership parameters

Controls how a `gene x k` loading matrix from ICA, NMF or DGRDL is
turned into module membership. Genes are kept where they sit in the tail
of a component's loading distribution, which means membership is not
exclusive: a gene loading strongly on three components belongs to three
modules. Genes in no tail belong to nothing, which is the background
category an argmax assignment cannot give you.

## Usage

``` r
params_module_membership(
  method = c("zscore", "fdr"),
  cutoff = 3,
  fdr = 0.05,
  tails = c("auto", "upper", "both")
)
```

## Arguments

- method:

  String. `"zscore"` standardises each component robustly (median and
  MAD) and keeps `abs(z) > cutoff`. `"fdr"` converts to two-sided
  p-values against a Normal null fitted the same way, Benjamini-Hochberg
  adjusts, and keeps `padj < fdr`. Defaults to `"zscore"`.

- cutoff:

  Float. Absolute z threshold for `method = "zscore"`. Defaults to
  `3.0`.

- fdr:

  Float. Adjusted p-value threshold for `method = "fdr"`. Defaults to
  `0.05`.

- tails:

  String. `"auto"` uses an upper-tail-only test when every loading is
  non-negative (the NMF case) and a two-sided one otherwise. `"upper"`
  and `"both"` force the choice. Defaults to `"auto"`.

## Value

A list with the parameters for usage in subsequent functions.

## References

Biton, et al., Cell Rep, 2014
