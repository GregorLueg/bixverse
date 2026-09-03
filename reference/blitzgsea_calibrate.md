# Calibrate the blitzGSEA null model for a signature

Draws random gene sets across a log-spaced grid of set sizes, fits gamma
tails to the resulting enrichment scores and smooths the fitted
parameters across sizes, see Lachmann, et al. No gene set library enters
here, so one calibration serves every library you score against that
signature. This is where nearly all of the runtime sits; scoring
afterwards costs one gamma tail evaluation per pathway.

The returned null is a plain list, so it survives
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and can be cached
against a signature.

## Usage

``` r
blitzgsea_calibrate(stats, blitz_params = params_blitzgsea())
```

## Arguments

- stats:

  Named numeric vector. The gene level statistic. Sorted internally, so
  the order you hand it in does not matter.

- blitz_params:

  List. The blitzGSEA parameters, see
  [`params_blitzgsea()`](https://gregorlueg.github.io/bixverse/reference/params_blitzgsea.md)
  wrapper function.

## Value

An object of class `BlitzGseaNull`, a list with the following elements:

- anchor_sizes - Numeric vector. The anchor set sizes, ascending.

- shape_pos - Numeric vector. Smoothed positive-tail gamma shape.

- scale_pos - Numeric vector. Smoothed positive-tail gamma scale.

- shape_neg - Numeric vector. Smoothed negative-tail gamma shape.

- scale_neg - Numeric vector. Smoothed negative-tail gamma scale.

- pos_ratio - Numeric vector. Smoothed fraction of positive null scores
  at each anchor.

- ks_pos - Float. Mean goodness-of-fit p-value for the positive tail.

- ks_neg - Float. Mean goodness-of-fit p-value for the negative tail.

- centred - Boolean. Whether the signature was centred.

- n_genes - Integer. Size of the signature it was calibrated on.

## References

Lachmann, et al., Bioinformatics, 2022
