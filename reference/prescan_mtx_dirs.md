# Prescan multiple mtx directories for a multi-load

Walks each input directory, reads the features file to build the
**intersection** of gene IDs across inputs (matched by the first column,
typically Ensembl gene IDs), decompresses any `.mtx.gz` files into a
temporary directory, and builds the file tasks expected by
[`load_multi_mtx()`](https://gregorlueg.github.io/bixverse/reference/load_multi_mtx.md).

Each input directory must contain the standard 10x trio: a `.mtx` (or
`.mtx.gz`) file, a barcodes file, and a features/genes file. File names
are matched by extension; if a directory contains multiple matching
files an error is raised.

## Usage

``` r
prescan_mtx_dirs(
  dirs,
  exp_ids,
  cells_as_rows = FALSE,
  has_hdr = FALSE,
  .verbose = TRUE
)
```

## Arguments

- dirs:

  Character vector of input directories. Length \>= 2.

- exp_ids:

  Character vector of experiment identifiers, one per directory. Must be
  unique.

- cells_as_rows:

  Boolean. Applied uniformly to all inputs. Defaults to `FALSE` (10x
  convention: genes are rows).

- has_hdr:

  Boolean. Whether the barcodes/features files have a header row.
  Applied uniformly. Defaults to `FALSE` (10x convention).

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with:

- universe - Character vector of gene IDs in the intersection, in the
  order they will appear in the final var table.

- universe_size - Length of the universe.

- file_tasks - List of per-input task lists for Rust and DuckDB.

- temp_files - Character vector of temp files created during
  decompression; the caller should
  [`unlink()`](https://rdrr.io/r/base/unlink.html) these after use.

## Examples

``` r
# two CellRanger style directories reduced to their shared gene universe
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
dirs <- c(tempfile("cr_a"), tempfile("cr_b"))
for (d in dirs) {
  dir.create(d, recursive = TRUE)
  write_cellranger_output(
    d, data$counts, data$obs, data$var,
    rows = "cells", format_type = "csv", .verbose = FALSE
  )
}
scan_res <- prescan_mtx_dirs(
  dirs = dirs,
  exp_ids = c("a", "b"),
  cells_as_rows = TRUE,
  has_hdr = TRUE,
  .verbose = FALSE
)
scan_res$universe_size
#> [1] 40

unlink(c(dirs, scan_res$temp_files), recursive = TRUE, force = TRUE)
```
