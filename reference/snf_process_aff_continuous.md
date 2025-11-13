# Helper function to process continuous data for SNF

Helper function to process continuous data for SNF

## Usage

``` r
snf_process_aff_continuous(data, k, mu, distance_metric, normalise)
```

## Arguments

- data:

  Numerical matrix. Is of type samples x features and is named.

- k:

  Integer. Number of neighbours to consider.

- mu:

  Float. Normalisation factor for the Gaussian kernel width.

- distance_metric:

  String. One of `c("euclidean", "manhattan", "canberra", "cosine")`.
  Which distance metric to use for the continuous calculations.

- normalise:

  Boolean. Shall continuous values be Z-scored.

## Value

The affinity matrix based on continuous values.
