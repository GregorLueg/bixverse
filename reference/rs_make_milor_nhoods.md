# Generate the neighbourhoods akin to the miloR approach

**\[experimental\]**

## Usage

``` r
rs_make_milor_nhoods(
  embd,
  knn_indices,
  sample_ids,
  n_samples,
  milor_params,
  seed,
  verbose
)
```

## Arguments

- embd:

  Numeric matrix. Represents the matrix used to generate the kNN graph
  and will be used to refine the neighbourhoods.

- knn_indices:

  Integer matrix. Each row represents a given cell and the columns the
  neighbours. (0-indexed!)

- sample_ids:

  Integer vector. 0-indexed(!) sample label per cell, in `0..n_samples`.
  One entry per row of `embd`.

- n_samples:

  Integer. Number of distinct samples.

- milor_params:

  Named list. Contains the parameters for running the miloR approach.

- seed:

  Integer. Seed for reproducibility.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following elements:

- index_cell - Integer. 0-indexed positions of the cells defining the
  neighbourhood.

- nhoods_i - Integer. 0-indexed positions of the cells in the
  neighbourhood.

- nhoods_j - Integer. To which neighbourhood the cell belongs.

- nhoods_x - Numeric. The x-value of the COO type matrix, i.e., defaults
  to `1.0`.

- nrows - Integer. Number of cells in the matrix

- ncols - Integer. Number of refined neighbourhoods.

- kth_distances - The k-th distances for spatial FDR calculations.

- sample_counts - Numeric matrix of neighbourhoods x samples. The cells
  of each sample found in each neighbourhood.

- nhood_overlap - Numeric. Cells each neighbourhood shares with all the
  others, the `"graph-overlap"` weighting for the spatial FDR.
