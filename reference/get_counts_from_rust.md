# Helper function to get counts from Rust

Helper function to get counts from Rust

## Usage

``` r
get_counts_from_rust(
  rust_con,
  assay,
  return_format,
  cell_indices,
  gene_indices,
  .verbose = TRUE
)
```

## Arguments

- rust_con:

  `SingeCellCountData` class. The connector to Rust.

- assay:

  String. One of `c("raw", "norm")`.

- return_format:

  String. One of `c("cell", "gene")`.

- cell_indices:

  Optional integer vector. The index positions of cells to return.

- gene_indices:

  Optional integer vector. The index positions of genes to return.

- .verbose:

  Boolean. Controls verbosity

## Value

Returns a list with:

- indptr - The index pointers of the compressed sparse format.

- indices - The indices of the data.

- data - The underlying data.

- no_cells - Number of cells.

- no_genes - Number of genes.
