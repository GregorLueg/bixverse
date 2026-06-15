# Get the index position for a gene

Returns the index for a given gene based on the internal gene mapping.
This is used for the single cell-related classes and methods.

## Usage

``` r
get_cell_indices(x, cell_ids, rust_index)

## S7 method for class <bixverse::MetaCells>
get_cell_indices(x, cell_ids, rust_index)

# S3 method for class 'ScMap'
get_cell_indices(x, cell_ids, rust_index)

## S7 method for class <bixverse::SingleCells>
get_cell_indices(x, cell_ids, rust_index)
```

## Arguments

- x:

  An object to get the gene index from.

- cell_ids:

  String vector. The cell ids to search for.

- rust_index:

  Bool. Shall rust-based indexing be returned.

## Value

The indices of the cells
