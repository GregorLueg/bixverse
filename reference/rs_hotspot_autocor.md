# Calculate gene spatial auto-correlations

This function implements the HotSpot auto-correlation functionality and
will return to what extent a given gene shows auto-correlation in the
generated kNN-graph from the embeddings. For details see DeTomaso, et
al.

## Usage

``` r
rs_hotspot_autocor(
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

- gene_idx - 0-based integer indicating the gene index.

- gaerys_c - Gaery's C calculation for the autocorrelation coefficient.

- z_score - Z-score of the auto-correlation.

- pval - P-value derived from the Z-score.

- fdr - False discovery rate based on the p-value.

## References

DeTomaso, et al., Cell Systems, 2021
