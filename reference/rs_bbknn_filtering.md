# Reduce BBKNN results to Top X neighbours

**\[experimental\]** Keeps the first `no_neighbours_to_keep` stored
entries of each CSR row.

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

A list with `indices` and `dist`, both numeric (double) matrices of
shape (n_cells, no_neighbours_to_keep). Positions without neighbours are
`NaN` in both.
