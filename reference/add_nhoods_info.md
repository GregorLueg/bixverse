# Add neighbourhood info on majority cell type

This function adds cell type composition information to the
`nhoods_info` slot within the `sc_miloR` object. For each neighbourhood,
it calculates the proportion of the majority cell type and identifies
which cell type is most abundant. This is useful for annotating
differential abundance results with the cellular composition of each
neighbourhood.

## Usage

``` r
add_nhoods_info(x, cell_info)

# S3 method for class 'sc_miloR'
add_nhoods_info(x, cell_info)
```

## Arguments

- x:

  `sc_miloR` object on which to tag on additional neighbourhood
  information.

- cell_info:

  Character vector. Represents the cell type annotations you wish to add
  to the different neighbourhoods. Must be the same length as the number
  of cells (rows) in the nhoods matrix.

## Value

Modified `sc_miloR` object with updated `nhoods_info` containing
`majority_celltype` and `majority_prop` columns.
