# Calculate batch LISI scores

Computes the Local Inverse Simpson's Index on batch labels using the kNN
graph. Measures the effective number of batches in each cell's
neighbourhood. Under perfect mixing LISI equals the number of batches,
under no mixing LISI equals 1.

## Usage

``` r
rs_batch_lisi(knn_mat, batch_vector)
```

## Arguments

- knn_mat:

  Integer matrix. The rows represent the cells and the columns the
  neighbour indices.

- batch_vector:

  Integer vector. The integers indicate to which batch a given cell
  belongs.

## Value

A list with the following items

- per_cell - Per-cell LISI scores

- mean_lisi - Mean LISI

- median_lisi - Median LISI
