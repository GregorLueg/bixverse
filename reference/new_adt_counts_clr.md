# Generates a new `ADTCounts` class

This function generates a new `ADTCounts` class which uses CLR
normalisation under the hood. You have the choice between Seurat-style
CLR (no negative values) and normal CLR (allows negative values).

## Usage

``` r
new_adt_counts_clr(
  raw_counts,
  cell_info,
  seurat_clr = FALSE,
  clean_clr_counts = TRUE,
  percentile = 0.01
)
```

## Arguments

- raw_counts:

  Numeric matrix. The raw ADT counts.

- cell_info:

  Named integer vector. Output of
  [`get_cell_info()`](https://gregorlueg.github.io/bixverse/reference/get_cell_info.md).
  Defines as elements the cell indices (R-based) and as names the
  barcodes.

- seurat_clr:

  Boolean. Shall a Seurat-style CLR be applied.

- clean_clr_counts:

  Boolean. Shall the per-protein 1st percentile be removed from the CLR
  normalised counts.

- percentile:

  Numeric. The percentile to remove to reduce background effects.

## Value

`ADTCounts` that contains the raw and normalised ADT counts.

## Examples

``` r
# CLR normalisation of synthetic ADT counts
adt <- generate_single_cell_test_data_adt()
cell_info <- stats::setNames(
  seq_len(nrow(adt$counts)),
  rownames(adt$counts)
)
new_adt_counts_clr(adt$counts, cell_info = cell_info)
#> ADTCounts
#>   Cells:     1000 
#>   Proteins:  15 
#>   Type:      CLR 
```
