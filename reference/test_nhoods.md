# Test neighbourhoods for differential abundance

Performs differential abundance testing on single-cell neighbourhoods
using edgeR's quasi-likelihood negative binomial framework. The function
fits a generalised linear model to neighbourhood cell counts, tests for
differential abundance between conditions, and applies spatial FDR
correction to account for overlapping neighbourhoods. This
implementation follows the approach described in Dann et al., using
graph-based neighbourhoods to identify regions of significant
compositional changes in single-cell data.

## Usage

``` r
test_nhoods(
  x,
  design,
  design_df,
  coef = NULL,
  norm_method = c("TMM", "RLE", "logMS"),
  min_mean = 0,
  robust = TRUE,
  fdr_weighting = c("k-distance", "graph-overlap", "none")
)

# S3 method for class 'sc_miloR'
test_nhoods(
  x,
  design,
  design_df,
  coef = NULL,
  norm_method = c("TMM", "RLE", "logMS"),
  min_mean = 0,
  robust = TRUE,
  fdr_weighting = c("k-distance", "graph-overlap", "none")
)
```

## Arguments

- x:

  `sc_miloR` object for which to run the differential abundance
  analysis.

- design:

  Formula for the experimental design

- design_df:

  data.frame. Contains the metadata to be used for the generation of the
  model matrix.

- coef:

  Optional string/integer. For more complex experimental designs, you
  can specify which coefficient to test. If NULL, tests the last
  coefficient in the design matrix (typically the main effect of
  interest).

- norm_method:

  String. Normalisation method to use. One of
  `c("TMM", "RLE", "logMS")`. Defaults to TMM (trimmed mean of
  M-values).

- min_mean:

  Numeric. Minimum mean count threshold for filtering neighbourhoods.
  Neighbourhoods with mean counts below this value are excluded.
  Defaults to 0 (no filtering).

- robust:

  Logical. If TRUE, uses robust estimation of the quasi-likelihood
  dispersion. Recommended for datasets with potential outliers. Defaults
  to TRUE.

- fdr_weighting:

  String. Spatial FDR weighting scheme. One of
  `c("k-distance", "graph-overlap", "none")`. k-distance uses the
  distance to the k-th nearest neighbour, graph-overlap uses
  neighbourhood overlap counts. Defaults to k-distance.

## Value

The `sc_miloR` object with added model and results from the differential
abundance analysis.

## References

Dann et al., 2022, Nat Biotechnol
