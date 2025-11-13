# Reduce BBKNN matrix to Top X neighbours

Reduce BBKNN matrix to Top X neighbours

## Usage

``` r
rs_bbknn_filtering(indptr, indices, no_neighbours_to_keep)
```

## Arguments

- indptr:

  Integer vector. The index pointers of the underlying data.

- indices:

  Integer vector. The indices of the nearest neighbours.

- no_neighbours_to_keep:

  Integer. Number of nearest neighbours to keep.

## Value

A numerical matrix with the Top X neighbours per row. If
`no_neighbours_to_keep` is larger than the number of neighbours in the
data, these positions will be `NA`.
