# Generate SuperCells.

**\[experimental\]** This function implements the approach from Bilous,
et al. to generate meta cells or called here SuperCells. You can provide
pre-computed kNN data (indices + distances) via `knn_data`, or an
embedding via `embd` from which the kNN graph will be generated. You
need to at least provide `knn_data` or `embd`. When `cells_to_use` is
supplied, the kNN graph is always regenerated on the subset and any
`knn_data` is ignored. Distances are required when the SuperCell
parameters request the kernel-weighted graph.

## Usage

``` r
rs_supercell(
  f_path,
  embd,
  cells_to_keep,
  cells_to_use,
  knn_data,
  supercell_params,
  target_size,
  seed,
  verbose
)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- embd:

  Optional numerical matrix. The embedding matrix (for example PCA
  embedding) used for the generation of the kNN graph. Required when
  `knn_data` is not provided, and required when using `cells_to_use`.

- cells_to_keep:

  Optional indices of the cells to keep, i.e., the cells used for the
  generation of the embedding.

- cells_to_use:

  Optional indices of cells to use for meta cell generation. Useful if
  you wish to generate meta cells in specific cell types. If this is
  provided, `embd` and `cells_to_keep` are required and the kNN graph
  will be regenerated on the subset.

- knn_data:

  Optional list. This contains pre-computed kNN data (including
  distances). The user has to ensure consistency! Ignored when
  `cells_to_use` is set.

- supercell_params:

  A list containing the SuperCell parameters.

- target_size:

  Numeric. Target library size for re-normalisation of the meta cells.
  Typically `1e4`.

- seed:

  Integer. For reproducibility purposes.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following elements:

- assignments - A list containing assignment information with elements:
  assignments (vector), metacells (list), unassigned (vector),
  n_metacells, n_cells, n_unassigned

- aggregated - A list with indptr, indices, raw_counts, norm_counts,
  nrow, ncol in sparse format.
