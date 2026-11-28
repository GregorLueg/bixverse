# Read the cell barcodes Cell Ranger ATAC called

Read the cell barcodes Cell Ranger ATAC called

## Usage

``` r
read_cellranger_barcodes(path)
```

## Arguments

- path:

  String. `singlecell.csv` (rows with `is__cell_barcode == 1`) or a
  filtered `*_bc_matrix.h5` (its `matrix/barcodes`).

## Value

Character vector of barcodes, in file order.
