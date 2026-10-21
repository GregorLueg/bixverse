# Helper function to generate the kNN graph defaults

Controls how a kNN is turned into the undirected, weighted graph that
the clustering runs on. `"snn"` is the shared nearest neighbour graph of
bluster/Seurat, `"umap"` the fuzzy simplicial set ScanPy builds (its
`obsp/connectivities`) and `"knn"` the plain kNN union. Each method only
reads its own parameters.

## Usage

``` r
params_knn_graph_defaults()
```

## Value

A named list with the following elements:

- knn_sym - String. How to symmetrise the kNN into a graph. Clustering
  resolutions do not carry over between the methods: the `"umap"` and
  `"knn"` graphs only hold kNN pairs, the sNN graph also links cells
  that share neighbours and prunes weak edges. The same resolution can
  give more or fewer communities depending on the data, so retune it
  when switching. One of `c("snn", "umap", "knn")`. Defaults to `"snn"`.

- snn_similarity - String. `"snn"` only. `"jaccard"` weighs an edge by
  the Jaccard index of the two (self-inclusive) neighbour sets. `"rank"`
  takes the shared neighbour with the lowest combined rank, so closer
  shared neighbours give heavier edges. Both are scaled to `[0, 1]`. One
  of `c("jaccard", "rank")`. Defaults to `"jaccard"`.

- pruning - Numeric. `"snn"` only. Edges below this weight are dropped.
  For the Jaccard weight a threshold `p` keeps pairs that share at least
  a fraction `2p / (1 + p)` of their (self-inclusive) neighbourhoods,
  rounded up to whole cells, so the same value means the same thing at
  any `k`. The default asks for about 12%, i.e. 2 of 16 shared cells at
  `k = 15`, in line with Seurat's `1/15`. Avoid a threshold that equals
  one of the discrete weights exactly, or float rounding decides whether
  the edge survives. Defaults to `0.0625`.

- rescue - Boolean. `"snn"` only. Reconnect the small components that
  pruning cuts off via their strongest sNN edge, so they do not turn up
  as tiny clusters. Defaults to `TRUE`.

- local_connectivity - Numeric. `"umap"` only. Number of nearest
  neighbours assumed to be fully connected, as in umap-learn. Defaults
  to `1.0`.

- mix_ratio - Numeric. `"umap"` only. `1` is the fuzzy union of the
  directed kNN memberships, `0` the fuzzy intersection (umap-learn's
  `set_op_mix_ratio`). Defaults to `1.0`.

- binary - Boolean. `"knn"` only. If `TRUE`, every edge weighs 1.
  Otherwise mutual neighbours weigh 2 and one-sided ones 1. Defaults to
  `FALSE`.

## References

Xu & Su, Bioinformatics, 2015; McInnes, et al., arXiv, 2018.
