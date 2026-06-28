# Run scDblFinder doublet detection on a set of cells

Run scDblFinder doublet detection on a set of cells

## Usage

``` r
.scdblfinder_run(
  object,
  cells_to_use,
  scdblfinder_params,
  return_features,
  streaming,
  seed,
  .verbose
)
```

## Arguments

- object:

  A `SingleCells` object.

- cells_to_use:

  Integer vector of 0-indexed cell indices.

- scdblfinder_params:

  List of scDblFinder parameters from
  [`params_scdblfinder()`](https://gregorlueg.github.io/bixverse/reference/params_scdblfinder.md).

- return_features:

  Logical. If `TRUE`, the classifier feature matrix is included in the
  result, with cells as rows and feature names as column names.

- streaming:

  Optional logical. Whether to stream the count data. If `NULL`,
  resolved automatically via
  [`auto_streaming()`](https://gregorlueg.github.io/bixverse/reference/auto_streaming.md).

- seed:

  Integer. Random seed.

- .verbose:

  Logical or integer. Verbosity level.

## Value

A `ScDblFinderRes` object with `cell_indices` set as an attribute.
