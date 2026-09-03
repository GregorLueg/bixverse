# Bixverse implementation of the blitzGSEA algorithm

Rust-based version of blitzGSEA, see Lachmann, et al. Instead of
permuting per pathway, it calibrates a gamma null once for the signature
and then reads each pathway's p-value straight off the fitted tail. That
makes the cost independent of how many pathways you throw at it, which
is where it pulls away from permutation methods on large libraries.

Hand back a `null_model` from
[`blitzgsea_calibrate()`](https://gregorlueg.github.io/bixverse/reference/blitzgsea_calibrate.md)
to reuse a calibration across libraries. Left as `NULL`, one is
calibrated on the fly.

## Usage

``` r
calc_blitzgsea(
  stats,
  pathways,
  blitz_params = params_blitzgsea(),
  null_model = NULL
)
```

## Arguments

- stats:

  Named numeric vector. The gene level statistic.

- pathways:

  List. A named list with each element containing the genes for this
  pathway.

- blitz_params:

  List. The blitzGSEA parameters, see
  [`params_blitzgsea()`](https://gregorlueg.github.io/bixverse/reference/params_blitzgsea.md)
  wrapper function. This function generates a list containing:

  - min_size - Integer. Minimum size for the gene sets.

  - max_size - Integer. Maximum size for the gene sets.

  - permutations - Integer. Random gene sets per anchor size.

  - anchors - Integer. Number of log-spaced anchor sizes.

  - symmetric - Boolean. Pool both tails into a single gamma.

  - centre - Boolean. Centre the signature before scoring.

  - ks_test - Boolean. Run the goodness-of-fit diagnostic.

  - seed - Float. Random seed.

- null_model:

  Optional `BlitzGseaNull` object from
  [`blitzgsea_calibrate()`](https://gregorlueg.github.io/bixverse/reference/blitzgsea_calibrate.md).
  Has to have been calibrated on the same signature and with the same
  `centre` setting. Defaults to `NULL`, in which case it is calibrated
  here.

## Value

A data.table with the results from the blitzGSEA with the following
columns:

- pathway_name - Character. The name of the pathway.

- es - Float. The enrichment score for this pathway.

- nes - Float. The normalised enrichment score, the signed normal
  quantile of the one-sided tail probability.

- pvals - Float. The two-sided p-value from the gamma approximation.

- fdr - Float. The Benjamini-Hochberg adjusted p-value.

- sidak - Float. The Sidak adjusted p-value.

- size - Integer. The size of the pathway after intersection.

- leading_edge - List of character vectors with the leading edge genes.

## References

Lachmann, et al., Bioinformatics, 2022
