# bixverse single cell ATAC class

A single cell ATAC experiment without RNA: a peak or genome tile store,
with the same on-disk layout as
[`SingleCells()`](https://gregorlueg.github.io/bixverse/reference/SingleCells.md).
The var table holds the features as `gene_id` (`chr:start-end`, 0-based
half-open) with `chrom`, `start` and `end`, so the usual store machinery
reads it unchanged. Create one with
[`load_atac_fragments()`](https://gregorlueg.github.io/bixverse/reference/load_atac_fragments.md)
or
[`load_atac_tenx()`](https://gregorlueg.github.io/bixverse/reference/load_atac_tenx.md),
reopen one with
[`load_existing()`](https://gregorlueg.github.io/bixverse/reference/load_existing.md).

`modality = "atac"` and the default `"rna"` both name this object's own
store and cache, so code written for the ATAC side of a
[`SingleCellsMultiModal()`](https://gregorlueg.github.io/bixverse/reference/SingleCellsMultiModal.md)
runs unchanged.

## Usage

``` r
SingleCellsAtac(dir_data)
```

## Arguments

- dir_data:

  String. The directory with the experimental files.

## Value

Returns the `SingleCellsAtac` class for further operations.

## Properties

- db_connection:

  R6 class with the DuckDB connection.

- count_connection:

  Rust connection to the peak or tile store.

- dir_data:

  Path to the directory with the data.

- sc_cache:

  The ATAC cache: LSI and other embeddings, kNN, graph.

- sc_map:

  Cell and feature mappings, the cell filter.

- dims:

  Dimensions of the store, cells x features.

- atac_info:

  `AtacInfo` with the store state and the LSI fit.

## Examples

``` r
# an empty handle; fill it with load_atac_fragments() or load_atac_tenx()
dir <- tempfile("bixverse_atac")
dir.create(dir)
object <- SingleCellsAtac(dir_data = dir)

unlink(dir, recursive = TRUE, force = TRUE)
```
