# Calculate gene to gene spatial correlations

**\[experimental\]** This function implements the HotSpot gene \<\> gene
local correlation functionality from HotSpot, see DeTomaso, et al.

## Usage

``` r
rs_hotspot_gene_cor(
  f_path_genes,
  f_path_cells,
  embd,
  knn_data,
  hotspot_params,
  cells_to_keep,
  genes_to_use,
  working_mem_gb,
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

- knn_data:

  Optional list. This contains pre-computed kNN data (including
  distances). The user has to ensure consistency! If provided, this will
  be used.

- hotspot_params:

  List. The HotSpot parameter list.

- cells_to_keep:

  Integer vector. 0-index vector indicating which cells to include in
  the analysis. Ensure that this is of same order/length as the
  embedding matrix.

- genes_to_use:

  Integer vector. 0-index vector indicating which genes to include.

- working_mem_gb:

  Numeric. Approximate working memory (GB) the streaming pair path may
  use for resident gene panels. Ignored when `streaming` is `FALSE`.
  Larger values mean fewer disk re-reads. Note this excludes the two
  dense N_genes x N_genes output matrices, which scale with
  `genes_to_use`.

- streaming:

  Boolean. Shall the data be streamed in chunks. Useful for large data
  sets.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

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
