# Write a `SingleCells` object to h5ad

Writes the counts, the DuckDB obs/var tables and whatever is cached in
memory (PCA, embeddings, neighbour graph) to a spec-compliant h5ad file,
so that the experiment can be handed to ScanPy or shared. The counts are
streamed cell batch by cell batch and never fully materialised in R.

Only the cells that are currently kept are written, i.e. the export
matches what `object[]` and `object[[]]` return. Counts are stored as
`float32`, which is the ScanPy convention. A failed write leaves no file
behind.

A `SingleCellsMultiModal` object inherits this method and exports its
RNA modality; the ADT and ATAC layers are not written.

## Usage

``` r
save_h5ad(
  object,
  h5_path,
  assay = c("raw", "norm"),
  genes = NULL,
  chunk_size = 10000L,
  overwrite = TRUE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- h5_path:

  File path to write the h5ad file to.

- assay:

  String. One of `c("raw", "norm")`. Which count assay to place in `X`.

- genes:

  Optional character vector of gene ids or integer vector of (1-based)
  gene indices to write, in the column order wanted. `NULL` writes every
  gene.

- chunk_size:

  Integer. Number of cells per streaming batch. Defaults to `10000L`.

- overwrite:

  Boolean. Shall an existing file be overwritten. Defaults to `TRUE`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

Returns the path to the written file, invisibly.
