# Helper function to generate kNN data with distances

Wrapper class that stores kNN data for subsequent usage in various other
functions, methods and helpers.

## Usage

``` r
new_sc_knn(knn_data, used_cells)
```

## Arguments

- knn_data:

  Named list with the following items:

  - indices - Integer matrix containing the indices of the nearest
    neighbours (0-indexed).

  - dist - Numerical matrix containing the distances to the nearest
    neighbours.

  - dist_metric - String. Distance metric used.

- used_cells:

  Character vector. The cells used to generate the kNN graph with the
  distances.

## Value

Generates the `sc_knn` class.
