# Getter the var table

Getter the var table

## Usage

``` r
get_sc_var(object, indices = NULL, cols = NULL, modality = c("rna", "adt"))
```

## Arguments

- object:

  `SingleCells`, `MetaCells`, `SingleCellsMultiModal` class.

- indices:

  Optional integer vector. The integer positions of the genes to return.

- cols:

  Optional string vector. The columns from the var table to return.

- modality:

  String. The modality to return. One of `c("rna", "adt")`.

## Value

The vars table
