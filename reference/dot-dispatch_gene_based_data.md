# Dispatch CSR-to-CSC generation based on streaming level

Internal helper to keep the streaming dispatch consistent across all
loaders. Validates the streaming level and routes to the appropriate
Rust method on the count connector.

## Usage

``` r
.dispatch_gene_based_data(
  rust_con,
  streaming,
  batch_size,
  max_genes_in_memory,
  cell_batch_size,
  .verbose
)
```

## Arguments

- rust_con:

  The Rust count connector.

- streaming:

  Integer. `0L` for in-memory, `1L` for light streaming, `2L` for heavy
  streaming with memory upper boundaries.

- batch_size:

  Integer. Batch size for light streaming.

- max_genes_in_memory:

  Integer. Genes held in memory at once for heavy streaming.

- cell_batch_size:

  Integer. Cell batch size for heavy streaming.

- .verbose:

  Boolean.

## Value

Invisible NULL. Side effect is the gene-based binary file.
