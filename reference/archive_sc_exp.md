# Archive the on-disk counts for cold storage

Compresses `counts_cells.bin` into `counts.bxa` with zstd. The
gene-based file is not archived, as it is a transpose of the cell-based
one. Normalised values are only kept where they cannot be recomputed
from the raw counts.
[`load_existing()`](https://gregorlueg.github.io/bixverse/reference/load_existing.md)
restores both binaries automatically when it finds only the archive.
Save the in-memory data first with
[`save_sc_exp_to_disk()`](https://gregorlueg.github.io/bixverse/reference/save_sc_exp_to_disk.md).

## Usage

``` r
archive_sc_exp(object, level = 3L, remove_bins = TRUE, .verbose = TRUE)
```

## Arguments

- object:

  `SingleCells` class.

- level:

  Integer. zstd compression level, between 1 and 22. Higher levels
  compress more and take longer; decompression speed barely changes.

- remove_bins:

  Boolean. Delete `counts_cells.bin` and `counts_genes.bin` once the
  archive is written.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with `n_cells`, `nnz`, `n_norm_stored` and `archive_bytes`,
invisibly.

## Examples

``` r
sc <- demo_single_cells()
dir <- sc@dir_data
save_sc_exp_to_disk(sc, type = "rds")
archive_sc_exp(sc, .verbose = FALSE)
sc <- load_existing(SingleCells(dir_data = dir), .verbose = FALSE)

unlink(dir, recursive = TRUE, force = TRUE)
```
