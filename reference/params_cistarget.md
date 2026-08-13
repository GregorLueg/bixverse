# Wrapper function to CisTarget parameters

`auc_threshold` and `max_rank` are two different cutoffs and are easy to
confuse. The first truncates the recovery curve used to score motif
enrichment, the second sets how deep into the ranking the background
curve (mean + 2 SD across all motifs) is built, and therefore where the
leading edge cuts. RcisTarget uses `maxRank = 5000` and `nMean = 100`;
the defaults here follow it. `auc_threshold` follows pySCENIC at 5%,
RcisTarget itself uses 3%.

## Usage

``` r
params_cistarget(
  auc_threshold = 0.05,
  nes_threshold = 3,
  max_rank = 5000L,
  n_mean = 100L,
  rcc_method = c("approx", "icistarget"),
  high_conf_cats = c("directAnnotation", "inferredBy_Orthology"),
  low_conf_cats = c("inferredBy_MotifSimilarity",
    "inferredBy_MotifSimilarity_n_Orthology")
)
```

## Arguments

- auc_threshold:

  Numeric between 0 and 1. Proportion of genes to use for AUC threshold
  calculation. Default is 0.05 (5% of genes).

- nes_threshold:

  Numeric. Normalised Enrichment Score threshold for significant motifs.
  Default is 3.0.

- max_rank:

  Integer. Depth of the recovery curves used to derive the background
  and the leading edge. Clamped to the number of genes in the ranking
  database. Default is 5000, the RcisTarget value.

- n_mean:

  Integer. Window for the rolling mean smoothing the background recovery
  curve. Only read when `rcc_method = "approx"`. Default is 100, the
  RcisTarget value.

- rcc_method:

  Character. Method for recovery curve calculation. Either "approx"
  (approximate, faster) or "icistarget" (exact, slower).

- high_conf_cats:

  Character vector. Annotation categories considered high confidence.
  Default includes direct annotations and orthology-based inferences.

- low_conf_cats:

  Character vector. Annotation categories considered lower confidence.
  Default includes motif similarity-based inferences.

## Value

A validated list of RcisTarget parameters.
