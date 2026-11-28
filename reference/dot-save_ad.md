# The shared body of `save_h5ad()` and `save_zarr()`

The shared body of
[`save_h5ad()`](https://gregorlueg.github.io/bixverse/reference/save_h5ad.md)
and
[`save_zarr()`](https://gregorlueg.github.io/bixverse/reference/save_zarr.md)

## Usage

``` r
.save_ad(object, path, zarr, assay, genes, chunk_size, overwrite, .verbose)
```

## Arguments

- object:

  `SingleCells` class.

- path:

  String. The h5ad file or Zarr directory to write.

- zarr:

  Boolean. Write an AnnData Zarr store instead of an h5ad file.

- assay:

  String. One of `c("raw", "norm")`.

- genes:

  Optional character vector of gene ids or integer vector of 1-based
  gene indices.

- chunk_size:

  Integer. Number of cells per streaming batch.

- overwrite:

  Boolean. Shall an existing output be overwritten.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

The written path, invisibly.
