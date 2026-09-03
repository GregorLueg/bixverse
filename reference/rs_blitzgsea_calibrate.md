# Calibrate the blitzGSEA gamma null for a signature

**\[experimental\]**

Draws random gene sets across a log-spaced grid of anchor sizes, fits
gamma tails to the resulting enrichment scores and smooths the fitted
parameters across sizes. Nothing about any gene set library enters here,
so one calibration serves every library scored against that signature.

## Usage

``` r
rs_blitzgsea_calibrate(stats, blitz_params)
```

## Arguments

- stats:

  Numeric vector. The gene level statistic. Needs to be sorted in
  descending nature.

- blitz_params:

  List. The blitzGSEA parameters, see
  [`params_blitzgsea()`](https://gregorlueg.github.io/bixverse/reference/params_blitzgsea.md).
  Recognised elements are `permutations`, `anchors`, `symmetric`,
  `centre`, `ks_test` and `seed`; anything else is ignored and any
  missing element takes its default.

## Value

List with the following elements

- anchor_sizes Numeric vector. The anchor set sizes, ascending.

- shape_pos Numeric vector. Smoothed positive-tail gamma shape.

- scale_pos Numeric vector. Smoothed positive-tail gamma scale.

- shape_neg Numeric vector. Smoothed negative-tail gamma shape.

- scale_neg Numeric vector. Smoothed negative-tail gamma scale.

- pos_ratio Numeric vector. Smoothed fraction of positive null scores at
  each anchor.

- ks_pos Float. Mean goodness-of-fit p-value for the positive tail.

- ks_neg Float. Mean goodness-of-fit p-value for the negative tail.

- centred Boolean. Whether the signature was centred.
