# Run Scrublet doublet detection on a set of cells

Run Scrublet doublet detection on a set of cells

## Usage

``` r
.scrublet_run(
  object,
  cells_to_use,
  scrublet_params,
  seed,
  streaming,
  return_combined_pca,
  return_pairs,
  .verbose
)
```

## Arguments

- object:

  A `SingleCells` object.

- cells_to_use:

  Integer vector of 0-indexed cell indices.

- scrublet_params:

  List of Scrublet parameters from
  [`params_scrublet()`](https://gregorlueg.github.io/bixverse/reference/params_scrublet.md).

- seed:

  Integer. Random seed.

- streaming:

  Optional logical. Whether to stream the count data. If `NULL`,
  resolved automatically via
  [`auto_streaming()`](https://gregorlueg.github.io/bixverse/reference/auto_streaming.md).

- return_combined_pca:

  Logical. If `TRUE`, PCA embeddings for observed cells and simulated
  doublets are included in the result.

- return_pairs:

  Logical. If `TRUE`, parent cell indices of simulated doublets are
  included in the result.

- .verbose:

  Logical or integer. Verbosity level.

## Value

A `ScrubletRes` object with `cell_indices` set as an attribute.
