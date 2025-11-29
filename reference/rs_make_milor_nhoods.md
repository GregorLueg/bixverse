# Generate the neighbourhoods akin to the miloR approach

Rust version of the

## Usage

``` r
rs_make_milor_nhoods(embd, knn_indices, milor_params, seed, verbose)
```

## Arguments

- embd:

  Numeric matrix. Represents the matrix used to generate the kNN graph
  and will be used to refine the neighbourhoods.

- knn_indices:

  Integer matrix. Each row represents a given cell and the columns the
  neighbours. (0-indexed!)

- milor_params:

  Named list. Contains the parameters for running the miloR approach.

- seed:

  Integer. Seed for reproducibility.

- verbose:

  Boolean. Controls verbosity of the function.

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
