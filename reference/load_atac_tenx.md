# Load an ATAC-only experiment from a 10x peak matrix h5

Load an ATAC-only experiment from a 10x peak matrix h5

## Usage

``` r
load_atac_tenx(
  object,
  h5_path,
  barcodes = NULL,
  max_mem_gb = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCellsAtac` class.

- h5_path:

  String. Cell Ranger ATAC's `filtered_peak_bc_matrix.h5`, or the Cell
  Ranger ARC h5 (its `Peaks` rows).

- barcodes:

  Optional character vector. Cells to load; `NULL` loads every barcode
  of the h5.

- max_mem_gb:

  Optional numeric. Memory budget of the feature-major conversion.

- .verbose:

  Boolean. Controls verbosity.

## Value

The loaded experiment, as
[`load_atac_fragments()`](https://gregorlueg.github.io/bixverse/reference/load_atac_fragments.md)
describes it; the counts are insertions.
