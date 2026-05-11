# Get the index position for a gene

Returns the index for a given gene based on the internal gene mapping.
This is used for the single cell-related classes and methods.

## Usage

``` r
get_gene_indices(x, gene_ids, rust_index)

# S3 method for class 'ScMap'
get_gene_indices(x, gene_ids, rust_index)
```

## Arguments

- x:

  An object to get the gene index from.

- gene_ids:

  String vector. The gene ids to search for.

- rust_index:

  Bool. Shall Rust-based indexing be returned.

## Value

The indices of the genes
