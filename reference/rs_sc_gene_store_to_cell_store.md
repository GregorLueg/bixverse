# Rebuilds the cell-major companion of a gene-major store

**\[experimental\]** Writes the `counts_cells.bin` twin of a
`counts_genes.bin` file. Memory is bounded by phasing over cells: each
phase holds one window of cells and re-reads the gene file to fill it,
so peak memory is the window rather than the matrix.

## Usage

``` r
rs_sc_gene_store_to_cell_store(
  f_path_in,
  f_path_out,
  cells_per_phase,
  gene_batch_size,
  verbose
)
```

## Arguments

- f_path_in:

  String. Path to the gene-major source file.

- f_path_out:

  String. Path of the cell-major file to write.

- cells_per_phase:

  Integer. Cells held in memory at once.

- gene_batch_size:

  Integer. Genes read per batch within a phase.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

String. The path that was written.
