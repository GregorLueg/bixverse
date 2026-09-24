# Wrapper function to generate SNF parameters

Wrapper function to generate SNF parameters

## Usage

``` r
params_snf(
  k = 20L,
  t = 20L,
  mu = 0.5,
  alpha = 1,
  normalise = TRUE,
  distance_metric = c("euclidean", "manhattan", "canberra", "cosine")
)
```

## Arguments

- k:

  Integer. Number of neighbours to consider. Defaults to `20L`.

- t:

  Integer. Number of iterations for the SNF algorithm. Defaults to
  `20L`.

- mu:

  Numeric. Normalisation factor for the Gaussian kernel width. Defaults
  to `0.5`.

- alpha:

  Numeric. Normalisation parameter controlling the fusion strength.
  Defaults to `1.0`.

- normalise:

  Boolean. Shall continuous values be Z-scored. Defaults to `TRUE`.

- distance_metric:

  String. Which distance metric to use for the continuous calculations.
  In case of pure categorical, Hamming will be used, for mixed data
  types Gower distance is used. One of
  `c("euclidean", "manhattan", "canberra", "cosine")`. Defaults to
  `"euclidean"`.

## Value

A named list with the following elements:

- k - Integer. Number of neighbours to consider. Defaults to `20L`.

- t - Integer. Number of iterations for the SNF algorithm. Defaults to
  `20L`.

- mu - Numeric. Normalisation factor for the Gaussian kernel width.
  Defaults to `0.5`.

- alpha - Numeric. Normalisation parameter controlling the fusion
  strength. Defaults to `1.0`.

- distance_metric - String. Which distance metric to use for the
  continuous calculations. In case of pure categorical, Hamming will be
  used, for mixed data types Gower distance is used. One of
  `c("euclidean", "manhattan", "canberra", "cosine")`. Defaults to
  `"euclidean"`.

- normalise - Boolean. Shall continuous values be Z-scored. Defaults to
  `TRUE`.
