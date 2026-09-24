# Calculate VISION pathway scores with auto-correlation (for meta cells)

**\[experimental\]** The function will take in a list of gene sets that
contains lists of `"pos"` and `"neg"` gene indices (0-indexed). You
don't have to provide the `"neg"`, but it can be useful to classify the
delta of two stats (EMT, Th1; Th2) etc. Additionally, it will take a
random gene list and calculate an auto-correlation score based on
Geary's C to identify pathways that show significant patterns on the kNN
graph generated on the provided embedding. This version works on
MetaCell counts which are stored in memory directly.

## Usage

``` r
rs_mc_vision_with_autocorrelation(
  sparse_data,
  embd,
  knn_data,
  gs_list,
  random_gs_list,
  vision_params,
  cluster_membership,
  verbose,
  seed
)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `cs_type`. Shape is (metacells, genes) and the data need to
  be the **normalised** counts.

- embd:

  Numerical matrix. The embedding matrix to use to generate the kNN
  graph. Needs to be of the same order/length as the meta cells in
  `sparse_data`.

- knn_data:

  Optional list. This contains pre-computed kNN data (`indices`, `dist`,
  `k`) and the `dist_metric` it was built with. The user has to ensure
  consistency! If provided, this will be used rather than a graph built
  from the parameter list.

- gs_list:

  Nested list. Each sublist contains the (0-indexed!) positive and
  negative gene indices of that specific gene set.

- random_gs_list:

  Double-nested list. The outer list represents the clusters of gene
  sets and the inner list represents the permutations within that
  cluster.

- vision_params:

  List. The kNN parameters, only read when no `knn_data` is provided.

- cluster_membership:

  Integer vector. 1-indexed(!) position in `random_gs_list` of the
  permuted cluster each gene set belongs to.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

- seed:

  Integer. Random seed for reproducibility.

## Value

A list with the following items:

- autocor_res - List with `auto_cor` (1 - Geary's C), `p_val` and `fdr`,
  one entry per gene set.

- vision_mat - A matrix of meta cells x vision scores per gene set.

## References

DeTomaso, et al., Nat. Commun., 2019
