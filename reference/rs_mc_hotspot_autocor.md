# Calculate gene spatial auto-correlations (for meta cells)

**\[experimental\]** This function implements the HotSpot
auto-correlation functionality and will return to what extent a given
gene shows auto-correlation in the kNN-graph over the meta cells. For
details see DeTomaso, et al. This version works on MetaCell counts which
are stored in memory directly. There is no streaming variant: streaming
bounds disk re-reads, which is not a problem an in-memory matrix has.

## Usage

``` r
rs_mc_hotspot_autocor(
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
  ensure consistency! If provided, this will be used and whether the
  distances are treated as squared is derived from `dist_metric` rather
  than from the parameter list.

- hotspot_params:

  List. The HotSpot parameter list. The kNN parameters are only read
  when no `knn_data` is provided.

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

- gene_idx - 0-based integer indicating the gene index.

- gaerys_c - Gaery's C calculation for the autocorrelation coefficient.

- z_score - Z-score of the auto-correlation.

- pval - P-value derived from the Z-score.

- fdr - False discovery rate based on the p-value.

## References

DeTomaso, et al., Cell Systems, 2021
