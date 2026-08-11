# Resolve cell names to kNN row positions

The kNN indices are positions within the cells the graph was built over,
not global cell indices, so cell names are resolved against
`used_cells`.

## Usage

``` r
.cells_to_knn_idx(cell_ids, used_cells, arg_name)
```

## Arguments

- cell_ids:

  Character vector. The cell names to resolve.

- used_cells:

  Character vector. The cells the kNN graph was built over.

- arg_name:

  String. Name of the argument, used in the error message.

## Value

Integer vector with the 0-indexed positions for Rust.
