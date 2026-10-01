# Build the node table of a Bonsai tree

Turns the 0-indexed parent array Rust hands back into a 1-indexed node
table, root parent `NA`. Leaves come first and carry the cell names.

## Usage

``` r
.bonsai_nodes(rs_res, cell_names, leaf_sizes)
```

## Arguments

- rs_res:

  List. Needs `parent`, `branch`, `x` and `y`.

- cell_names:

  Character vector. One name per leaf, in leaf order.

- leaf_sizes:

  Integer. Cells behind each leaf, same order.

## Value

data.table with `node`, `parent`, `branch`, `is_leaf`, `cell_id`,
`n_cells`, `x` and `y`.
