# The cells of an ATAC-only experiment, from exactly one source

The cells of an ATAC-only experiment, from exactly one source

## Usage

``` r
.resolve_atac_barcodes(
  fragments,
  barcodes,
  cellranger,
  min_fragments,
  exclude_chroms
)
```

## Arguments

- fragments:

  String. The fragment file.

- barcodes, cellranger, min_fragments:

  The sources, exactly one of them not `NULL`.

- exclude_chroms:

  Character vector. Chromosomes the `min_fragments` count skips.

## Value

Character vector of cell barcodes.
