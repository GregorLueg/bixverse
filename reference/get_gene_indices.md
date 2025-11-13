# Get the index position for a gene

Get the index position for a gene

Get the gene indices from a `single_cell_exp`.

## Usage

``` r
get_gene_indices(x, gene_ids, rust_index)

# S3 method for class 'sc_mapper'
get_gene_indices(x, gene_ids, rust_index)
```

## Arguments

- x:

  An object to get the gene index from.

- gene_ids:

  String vector. The gene ids to search for.

- rust_index:

  Bool. Shall rust-based indexing be returned.
