# Get the gene names based on the gene idx

Get the gene names based on the gene idx

Get the highly variable gene indices from a `single_cell_exp`.

## Usage

``` r
get_gene_names_from_idx(x, gene_idx, rust_based = TRUE)

# S3 method for class 'sc_mapper'
get_gene_names_from_idx(x, gene_idx, rust_based = TRUE)
```

## Arguments

- x:

  An object to get the gene names from.

- gene_idx:

  Integer. The original gene indices for which to return the gene names.

- rust_based:

  Boolean. Is it Rust-based, i.e., 0-index or R-based, i.e., 1-indexed.
