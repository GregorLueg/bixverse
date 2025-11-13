# Bixverse implementation of GSVA

Implementation of the bixverse version of the gene set variation
analysis (GSVA).

## Usage

``` r
calc_gsva(
  exp,
  pathways,
  gaussian = TRUE,
  gsva_params = params_gsva(),
  .verbose = FALSE
)
```

## Arguments

- exp:

  Numerical matrix. Rows represents the features, columns the
  features/genes.

- pathways:

  List. A named list with each element containing the genes for this
  pathway.

- gaussian:

  Boolean. If set to `TRUE` the Gaussian kernel will be used, if `FALSE`
  the Poisson will be used.

- gsva_params:

  List. The GSVA parameters, see [`params_gsva()`](params_gsva.md)
  wrapper function. This function generates a list containing:

  - tau - Float. The tau parameter of the algorithm. Large values will
    emphasise the tails more. Defaults to `1.0`.

  - min_size - Integer. Minimum size for the gene sets.

  - max_size - Integer. Maximum size for the gene sets.

  - max_diff - Boolean. Influences the scoring. If `TRUE` the difference
    will be used; if `FALSE`, the largest absolute value.

  - abs_rank - Boolean. If `TRUE` = pos-neg, `FALSE` = pos+neg for the
    internal calculations.

- .verbose:

  Boolean. Controls verbosity.

## Value

A matrix of shape pathways (that passed the thresholds) x samples.
