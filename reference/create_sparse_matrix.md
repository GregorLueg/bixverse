# Helper function transform Rust counts into sparse matrices

Helper function transform Rust counts into sparse matrices

## Usage

``` r
create_sparse_matrix(count_data, return_format)
```

## Arguments

- count_data:

  A list. Output of [`get_counts_from_rust()`](get_counts_from_rust.md).

- return_format:

  String. One of `c("cell", "gene")`.

## Value

The sparse matrix in CSR or CSC format, pending the choice in
`return_format`
