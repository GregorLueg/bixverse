# Generate meta cells (hdWGCNA method)

**\[experimental\]** This function implements the approach from
Morabito, et al. to generate meta cells. You can provide an already
pre-computed kNN matrix or an embedding to regenerate the kNN matrix
with specified parameters in the meta_cell_params. If `knn_mat` is
provided, this one will be used. You need to at least provide `knn_mat`
or `embd`!

## Usage

``` r
rs_get_metacells_bootstrapped(
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

  Optional integer vector. Original cell indices (0-indexed!) of the
  rows of `knn_mat` / `embd`, in row order. If `NULL`, the rows are
  assumed to map one-to-one onto the count file.

- cells_to_use:

  Optional integer vector. Original cell indices (0-indexed!) to
  restrict the meta cell generation to, e.g. specific cell types. Needs
  `cells_to_keep` and `embd`; the kNN graph is then regenerated on the
  subset. Cells not in `cells_to_keep` are dropped silently.

- meta_cell_params:

  A list containing the meta cell parameters.

- target_size:

  Numeric. Target library size for re-normalisation of the meta cells.
  Typically `1e4`.

- seed:

  Integer. For reproducibility.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following elements:

- assignments - A list with `assignments` (list with one integer vector
  of 1-indexed meta cell ids per cell, as meta cells can overlap),
  `metacells` (list of 1-indexed original cell indices per meta cell),
  `unassigned` (1-indexed cells in no meta cell), `n_metacells`,
  `n_cells` and `n_unassigned`.

- aggregated - A CSR list (meta cells x genes) with indptr, indices,
  raw_counts, norm_counts, nrow and ncol.
