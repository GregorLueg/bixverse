# Wrapper function to generate blitzGSEA parameters

Wrapper function to generate blitzGSEA parameters

## Usage

``` r
params_blitzgsea(
  min_size = 5L,
  max_size = 500L,
  permutations = 2000L,
  anchors = 40L,
  symmetric = FALSE,
  centre = TRUE,
  ks_test = TRUE,
  seed = 42
)
```

## Arguments

- min_size:

  Integer. Minimum number of genes per gene set.

- max_size:

  Integer. Maximum number of genes per gene set.

- permutations:

  Integer. Random gene sets drawn per anchor size during calibration.
  Defaults to `2000L`. Below `1000L` the two tails are pooled into a
  single gamma regardless of `symmetric`.

- anchors:

  Integer. Number of log-spaced anchor sizes requested. Sizes that
  collide after rounding are collapsed, so the realised grid is usually
  a little smaller. Defaults to `40L`.

- symmetric:

  Boolean. Pool both tails into one gamma instead of fitting them
  separately. Defaults to `FALSE`.

- centre:

  Boolean. Centre the signature on its mean before scoring. The
  enrichment score is not invariant to an offset, so the calibration and
  the scoring have to agree on this. Defaults to `TRUE`.

- ks_test:

  Boolean. Run the Kolmogorov-Smirnov goodness-of-fit diagnostic at
  every anchor. Costs a sort per anchor. Defaults to `TRUE`.

- seed:

  Float. Random seed for the calibration. Defaults to `42`.

## Value

List with parameters for usage in subsequent function.

## References

Lachmann, et al., Bioinformatics, 2022
