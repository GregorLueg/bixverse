# Restrict ADT counts to the kept cells

ADT counts are keyed by barcode and frozen when they are added, so they
can hold cells a later filter dropped. Returns the rows of the kept
cells present in the ADT matrix, in kept cell order.

## Usage

``` r
.adt_kept_rows(object, mat)
```

## Arguments

- object:

  `SingleCellsMultiModal` class.

- mat:

  Numeric matrix. Cells x features, cell ids as row names.

## Value

`mat` restricted to (and ordered by) the kept cells.
