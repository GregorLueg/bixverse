# bixverse meta cell class (deprecated)

**\[deprecated\]**

This constructor has been renamed to
[`MetaCells()`](https://gregorlueg.github.io/bixverse/reference/MetaCells.md).

## Usage

``` r
meta_cells(meta_cell_data, var_data, meta_cell_method)
```

## Arguments

- meta_cell_data:

  Named list. Output of meta-cell generation Rust functions. Will
  contain the aggregated raw and normalised counts, plus additional
  information on the origin of the meta cells.

- var_data:

  data.table with the variable/feature informations.

- meta_cell_method:

  String describing the origin of the metacell.

## Value

Returns a
[`MetaCells()`](https://gregorlueg.github.io/bixverse/reference/MetaCells.md)
object.
