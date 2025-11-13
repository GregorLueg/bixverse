# Calculate the SNF affinity matrix for continuous values

Calculate the SNF affinity matrix for continuous values

## Usage

``` r
rs_snf_affinity_continuous(data, distance_type, k, mu, normalise)
```

## Arguments

- data:

  Numerical matrix. Needs to be oriented features x samples!

- distance_type:

  String. One of `c("euclidean", "manhattan", "canberra", "cosine")`.
  Which distance metric to use here.

- k:

  Integer. Number of neighbours to consider.

- mu:

  Float. Normalisation factor for the Gaussian kernel width.

- normalise:

  Boolean. Shall continuous values be Z-scored.

## Value

The affinity matrix based on continuous values.
