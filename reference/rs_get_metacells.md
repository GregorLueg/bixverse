# Generate meta cells (hdWGCNA method)

This function implements the approach from Morabito, et al. to generate
meta cells. You can provide an already pre-computed kNN matrix or an
embedding to regenerate the kNN matrix with specified parameters in the
meta_cell_params. If `knn_mat` is provided, this one will be used. You
need to at least provide `knn_mat` or `embd`!

## Usage

``` r
rs_get_metacells(
  f_path,
  knn_mat,
  embd,
  cells_to_keep,
  cells_to_use,
  meta_cell_params,
  target_size,
  seed,
  verbose
)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- knn_mat:

  Optional integer matrix. The kNN matrix you wish to use for the
  generation of the meta cells. This function expects 0-indices!

- embd:

  Optional numerical matrix. The embedding matrix (for example PCA
  embedding) you wish to use for the generation of the kNN graph that is
  used subsequently for aggregation of the meta cells.

- cells_to_keep:

  Optional indices of the cells to keep, i.e., the cells used for the
  generation of the embedding.

- cells_to_use:

  Optional indices of cells to use for meta cell generation. Useful if
  you wish to generate meta cells in specific cell types. If this is
  provided, the kNN graph will be regenerated.

- meta_cell_params:

  A list containing the meta cell parameters.

- target_size:

  Numeric. Target library size for re-normalisation of the meta cells.
  Typically `1e4`.

- seed:

  Integer. For reproducibility purposes.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the following elements:

- assignments - A list containing assignment information with elements:
  assignments (vector), metacells (list), unassigned (vector),
  n_metacells, n_cells, n_unassigned

- aggregated - A list with indptr, indices, raw_counts, norm_counts,
  nrow, ncol in sparse format.
