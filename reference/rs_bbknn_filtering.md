# Reduce BBKNN results to the k nearest neighbours

**\[experimental\]** Keeps the batch-balanced `no_neighbours_to_keep`
nearest neighbours of each row of the BBKNN distance CSR: the nearest
`ceil(no_neighbours_to_keep / n_batches)` per batch, with slots a batch
cannot fill going to the nearest remaining neighbours. Rows come back
sorted by distance (ties by index). Rows with fewer stored entries are
padded by repeating their last neighbour; an empty row points at the
cell itself at distance zero.

## Usage

``` r
rs_bbknn_filtering(indptr, indices, data, batch_labels, no_neighbours_to_keep)
```

## Arguments

- indptr:

  Integer vector. The index pointers of the underlying data.

- indices:

  Integer vector. The indices of the nearest neighbours.

- data:

  Numeric vector. The distances to the nearest neighbours.

- batch_labels:

  Integer vector. The batch of each cell, as passed to
  [`rs_bbknn()`](https://gregorlueg.github.io/bixverse/reference/rs_bbknn.md).

- no_neighbours_to_keep:

  Integer. Number of nearest neighbours to keep.

## Value

A list with `indices` (integer, 0-indexed) and `dist` (numeric) matrices
of shape (n_cells, no_neighbours_to_keep), nearest first.
