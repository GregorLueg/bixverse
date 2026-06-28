# Identifies genes to include into a SCENIC analysis

**\[experimental\]**

## Usage

``` r
rs_scenic_gene_filter(f_path_genes, cell_indices, scenic_params, verbose)
```

## Arguments

- f_path_genes:

  Path to the `counts_genes.bin` file.

- cell_indices:

  Integer vector. 0-indexed(!) positions of cells to include in the
  analysis

- scenic_params:

  Named list. Contains all of the parameters need for SCENIC.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

The 0-indexed positions of the genes to include in the scenic analysis.
