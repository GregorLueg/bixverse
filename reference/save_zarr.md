# Write a `SingleCells` object to an AnnData Zarr store

Writes what
[`save_h5ad()`](https://gregorlueg.github.io/bixverse/reference/save_h5ad.md)
writes, as an AnnData Zarr (v3) directory that `anndata.read_zarr()`
opens. A failed write leaves no store behind.

Not available on Windows, where it errors: the Zarr writer pulls in
Blosc, whose bundled pthread shim clashes with the one Rust links there.
Use
[`save_h5ad()`](https://gregorlueg.github.io/bixverse/reference/save_h5ad.md)
instead. Reading Zarr stores works on Windows.

## Usage

``` r
save_zarr(
  object,
  zarr_path,
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

- zarr_path:

  Directory to write the Zarr store to.

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

Returns the path to the written store, invisibly.
