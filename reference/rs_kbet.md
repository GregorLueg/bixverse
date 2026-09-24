# Calculate kBET type scores

**\[experimental\]** The function takes in a kNN matrix and a batch
vector indicating which cell belongs to which batch. For the
neighbourhood of each cell, a chi-square test checks whether the batch
proportions differ from the overall batch proportions (with Yates'
correction for two batches). Good mixing means few cells with
significant differences; bad mixing means many.

## Usage

``` r
rs_kbet(knn_mat, batch_vector, verbose)
```

## Arguments

- knn_mat:

  Integer matrix. The rows represent the cells and the columns the
  neighbour indices (0-indexed!).

- batch_vector:

  Integer vector. The batch per cell. The codes need not be 0-based or
  contiguous.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the following items

- pval - Per-cell p-values from the chi-square test.

- chi_square_stats - Per-cell chi-square statistics.

- mean_chi_square - The mean chi-square value.

- median_chi_square - The median chi-square value.
