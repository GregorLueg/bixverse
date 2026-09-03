# Calculate gene to gene spatial correlations (for meta cells)

**\[experimental\]** This function implements the HotSpot gene \<\> gene
local correlation functionality from HotSpot, see DeTomaso, et al. This
version works on MetaCell counts which are stored in memory directly.

Three dense metacells x genes blocks are live at once, so keep
`genes_to_use` to the panel actually of interest rather than the whole
transcriptome.

## Usage

``` r
rs_mc_hotspot_gene_cor(
  sparse_data,
  embd,
  knn_data,
  hotspot_params,
  cells_to_keep,
  genes_to_use,
  verbose,
  seed
)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `format`. Shape is (metacells, genes) and the data are the
  raw counts.

- embd:

  Numerical matrix. The embedding matrix from which to generate the kNN
  graph.

- knn_data:

  Optional list. This contains pre-computed kNN data (including
  distances) and the `dist_metric` it was built with. The user has to
  ensure consistency! If provided, this will be used rather than a graph
  built from the parameter list.

- hotspot_params:

  List. The HotSpot parameter list. The kNN parameters are only read
  when no `knn_data` is provided; `normalise` is unused on this path.

- cells_to_keep:

  Integer vector. 0-index vector indicating which meta cells to include
  in the analysis. Ensure that this is of same order/length as the
  embedding matrix.

- genes_to_use:

  Integer vector. 0-index vector indicating which genes to include.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

- seed:

  Integer. Random seed for reproducibility.

## Value

A list with the following elements.

- cor - The gene x gene local correlation matrix.

- z - The Z-scores of these local correlations.

## References

DeTomaso, et al., Cell Systems, 2021
