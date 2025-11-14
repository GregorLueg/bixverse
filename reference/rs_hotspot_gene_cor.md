# Calculate gene\<\>gene spatial correlations

This function implements the HotSpot gene \<\> gene local correlation
functionality from HotSpot, see DeTomaso, et al.

## Usage

``` r
rs_hotspot_gene_cor(
  f_path_genes,
  f_path_cells,
  embd,
  hotspot_params,
  cells_to_keep,
  genes_to_use,
  streaming,
  verbose,
  seed
)
```

## Arguments

- f_path_genes:

  Path to the `counts_genes.bin` file.

- f_path_cells:

  Path to the `counts_cells.bin` file.

- embd:

  Numerical matrix. The embedding matrix from which to generate the kNN
  graph.

- hotspot_params:

  List. The HotSpot parameter list.

- cells_to_keep:

  Integer vector. 0-index vector indicating which cells to include in
  the analysis. Ensure that this is of same order/length as the
  embedding matrix.

- genes_to_use:

  Integer vector. 0-index vector indicating which genes to include.

- streaming:

  Boolean. Shall the data be streamed in chunks. Useful for large data
  sets.

- verbose:

  Boolean. Controls verbosity of the function.

- seed:

  Integer. Random seed for reproducibility.

## Value

A list with the following elements.

- cor - A matrix of the N x N genes_to_use length with the auto-
  correlation coefficients.

- z - A matrix of N x N genes_to_use length with the Z-scores of the
  local correlations between two genes.

## References

DeTomaso, et al., Cell Systems, 2021
