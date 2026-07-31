# Map 1-based ScType indices onto cell type names

Map 1-based ScType indices onto cell type names

## Usage

``` r
.sc_type_labels(idx, cell_types)
```

## Arguments

- idx:

  Integer vector. 1-based indices, `0L` denoting Unknown.

- cell_types:

  Character vector. The cell type names.

## Value

A character vector with `NA_character_` for the Unknowns.
