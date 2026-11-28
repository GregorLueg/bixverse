# Read counts from the ATAC store of a multi-modal object

The ATAC twin of the RNA path of
[`get_sc_counts()`](https://gregorlueg.github.io/bixverse/reference/get_sc_counts.md):
same cell filter, same matrix shapes, the features named by their
`feature_id`.

## Usage

``` r
.get_atac_counts(
  object,
  assay,
  return_format,
  cell_indices,
  gene_indices,
  use_cells_to_keep,
  .verbose
)
```

## Arguments

- object:

  `SingleCellsMultiModal` class with ATAC.

- assay:

  String. One of `c("raw", "norm")`.

- return_format:

  String. One of `c("cell", "gene")`.

- cell_indices:

  Optional integer vector of 1-based cell indices.

- gene_indices:

  Optional integer vector of 1-based feature indices.

- use_cells_to_keep:

  Boolean. Restrict to the kept cells.

- .verbose:

  Boolean. Controls verbosity.

## Value

A sparse matrix, as
[`get_sc_counts()`](https://gregorlueg.github.io/bixverse/reference/get_sc_counts.md).
