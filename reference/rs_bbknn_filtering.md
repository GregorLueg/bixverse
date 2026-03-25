# Reduce BBKNN results to Top X neighbours

Reduce BBKNN results to Top X neighbours

## Usage

``` r
rs_bbknn_filtering(indptr, indices, data, no_neighbours_to_keep)
```

## Arguments

- indptr:

  Integer vector. The index pointers of the underlying data.

- indices:

  Integer vector. The indices of the nearest neighbours.

- data:

  Numeric vector. The distances to the nearest neighbours.

- no_neighbours_to_keep:

  Integer. Number of nearest neighbours to keep.

## Value

A list with `indices` (integer matrix) and `dist` (numeric matrix), each
with shape (n_cells, no_neighbours_to_keep). Positions without
neighbours are filled with -1 (indices) or NaN (distances).
