# Generate SEACells

This function implements the SEACells algorithm for generating meta
cells from Persad et al. An embedding matrix must be provided which is
used to construct the kNN graph and kernel matrix for the SEACells
algorithm. This version is highly memory and speed-optimised and will
truncate small values during matrix operations which can affect
convergence.

## Usage

``` r
rs_get_seacells(
  f_path,
  embd,
  cells_to_keep,
  cells_to_use,
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

  Optional indices of the cells to keep, i.e., the cells used for the
  generation of the embedding.

- cells_to_use:

  Optional indices of cells to use for meta cell generation. Useful if
  you wish to generate meta cells in specific cell types.

- seacells_params:

  A list containing the SEACells parameters.

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

- rss - Vector of RSS values from each iteration.

- archetypes - Vector of cell indices selected as archetypes.

## References

Persad, et al., Nat. Biotechnol., 2023.
