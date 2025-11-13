# Wrapper function to CisTarget parameters

Wrapper function to CisTarget parameters

## Usage

``` r
params_cistarget(
  auc_threshold = 0.05,
  nes_threshold = 3,
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
