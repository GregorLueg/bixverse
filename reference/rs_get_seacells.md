# Generate SEACells

**\[experimental\]** This function implements the SEACells algorithm for
generating meta cells from Persad et al. An embedding matrix must be
provided which is used to construct the kNN graph and kernel matrix for
the SEACells algorithm. This version is highly memory and
speed-optimised and will truncate small values during matrix operations
which can affect convergence.

## Usage

``` r
rs_get_seacells(
  f_path,
  embd,
  cells_to_keep,
  cells_to_use,
  knn_data,
  seacells_params,
  target_size,
  seed,
  verbose
)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- embd:

  Numerical matrix. The embedding matrix (for example PCA embedding)
  used for the generation of the kNN graph and kernel matrix.

- cells_to_keep:

  Optional integer vector. Original cell indices (0-indexed!) of the
  rows of `embd`, in row order. If `NULL`, the rows are assumed to map
  one-to-one onto the count file.

- cells_to_use:

  Optional integer vector. Original cell indices (0-indexed!) to
  restrict the meta cell generation to, e.g. specific cell types. The
  kNN graph is then regenerated on the subset. Cells not in
  `cells_to_keep` are dropped silently.

- knn_data:

  Optional list. This contains pre-computed kNN data (including
  distances). The user has to ensure consistency! Ignored when
  `cells_to_use` is set.

- seacells_params:

  A list containing the SEACells parameters.

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
  Empty archetypes are dropped and the ids renumbered.

- aggregated - A CSR list (meta cells x genes) with indptr, indices,
  raw_counts, norm_counts, nrow and ncol.

- rss - Numerical vector of RSS values from each iteration.

- archetypes - Integer vector with the original cell indices
  (0-indexed!) of the archetypes of the retained meta cells.

## References

Persad, et al., Nat. Biotechnol., 2023.
