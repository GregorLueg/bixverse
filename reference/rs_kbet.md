# Calculate kBET type scores

The function takes in a kNN matrix and a batch vector indicating which
cell belongs to which batch. The function will check for the
neighbourhood of each cell if the proportion of represented batches are
different from the overall batch proportions. Good mixing of batches
would mean very cells have significant differences; bad mixing a lot of
the batches have bad mixing.

## Usage

``` r
rs_kbet(knn_mat, batch_vector)
```

## Arguments

- knn_mat:

  Integer matrix. The rows represent the cells and the columns the
  neighbour indices.

- batch_vector:

  Integer vector. The integers indicate to which batch a given cell
  belongs.

## Value

A vector of p-values based on the ChiSquare statistic per cell.
