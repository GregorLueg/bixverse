# Calculate the pairwise column distance in a matrix

This function allows to calculate pairwise between all columns the
specified distance metric.

## Usage

``` r
rs_dist(x, distance_type)
```

## Arguments

- x:

  Numerical matrix. The matrix for which to calculate the pairwise
  column distances.

- distance_type:

  String. One of `c("euclidean", "manhattan", "canberra", "cosine")`.

## Value

The calculated distance matrix
