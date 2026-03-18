# Get the index position for a gene

Get the index position for a gene

Set the gene mapping for a `SingleCells` class.

## Usage

``` r
get_cell_indices(x, cell_ids, rust_index)

# S3 method for class 'ScMap'
get_cell_indices(x, cell_ids, rust_index)
```

## Arguments

- x:

  An object to get the gene index from.

- cell_ids:

  String vector. The cell ids to search for.

- rust_index:

  Bool. Shall rust-based indexing be returned.
