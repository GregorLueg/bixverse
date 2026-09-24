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

  Optional integer vector. Original cell indices (0-indexed!) of the
  rows of `embd` / the kNN data, in row order. If `NULL`, the rows are
  assumed to map one-to-one onto the count file.

- cells_to_use:

  Optional integer vector. Original cell indices (0-indexed!) to
  restrict the meta cell generation to, e.g. specific cell types. If
  this is provided, `embd` and `cells_to_keep` are required and the kNN
  graph will be regenerated on the subset. Cells not in `cells_to_keep`
  are dropped silently.

- knn_data:

  Optional list. This contains pre-computed kNN data (including
  distances). The user has to ensure consistency! Ignored when
  `cells_to_use` is set.

- supercell_params:

  A list containing the SuperCell parameters. The number of meta cells
  is `ceiling(n_cells / graining_factor)`.

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

- assignments - A list with `assignments` (integer vector with the
  1-indexed meta cell id per original cell, `-1` if unassigned),
  `metacells` (list of 1-indexed original cell indices per meta cell),
  `unassigned` (1-indexed), `n_metacells`, `n_cells` and `n_unassigned`.

- aggregated - A CSR list (meta cells x genes) with indptr, indices,
  raw_counts, norm_counts, nrow and ncol.
