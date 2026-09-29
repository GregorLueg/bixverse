# Generate the gene-based binary from the cell-based one

Internal helper so every loader runs the CSR to CSC conversion, and
warns about the retired streaming arguments, the same way.

## Usage

``` r
.dispatch_gene_based_data(
  rust_con,
  csc_mem_gb,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose
)
```

## Arguments

- rust_con:

  The Rust count connector.

- csc_mem_gb:

  Optional numeric. Memory in GB for the conversion buffers. `NULL`
  converts in one pass.

- streaming, batch_size, max_genes_in_memory, cell_batch_size:

  Retired arguments, forwarded from the caller only to warn if they were
  supplied.

- .verbose:

  Boolean.

## Value

Invisible NULL. Side effect is the gene-based binary file.
