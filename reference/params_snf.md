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

  Integer. Number of neighbours to consider.

- t:

  Integer. Number of iterations for the SNF algorithm.

- mu:

  Float. Normalisation factor for the Gaussian kernel width.

- alpha:

  Float. Normalisation parameter controlling the fusion strength.

- normalise:

  Boolean. Shall continuous values be Z-scored.

- distance_metric:

  String. One of `c("euclidean", "manhattan", "canberra", "cosine")`.
  Which distance metric to use for the continuous calculations. In case
  of pure categorical, Hamming will be used, for mixed data types Gower
  distance is used.

## Value

List with parameters for usage in subsequent function.
