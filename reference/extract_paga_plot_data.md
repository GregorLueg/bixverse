# Extract the PAGA graph positioned on an embedding

Puts every cluster of a
[`run_paga_sc()`](https://gregorlueg.github.io/bixverse/reference/run_paga_sc.md)
result at the centroid of its cells in an embedding and returns the
abstracted graph as node and edge tables. A node-link plot drawn from
these sits where the reader already knows the biology is, which a free
layout of the abstracted graph does not.

The abstracted graph is close to complete on real data, so `threshold`
is doing real work. Drop it to zero and you get a hairball. `tree_only`
is the shortcut to the backbone.

## Usage

``` r
extract_paga_plot_data(
  object,
  paga_res,
  embedding = "umap",
  cluster_col = NULL,
  threshold = 0.01,
  tree_only = FALSE,
  centroid = c("median", "mean"),
  node_stat_col = NULL,
  ...
)
```

## Arguments

- object:

  A single cell class.

- paga_res:

  `PagaRes` class. The output of
  [`run_paga_sc()`](https://gregorlueg.github.io/bixverse/reference/run_paga_sc.md),
  run on this object.

- embedding:

  String. Name of the embedding to position the nodes in.

- cluster_col:

  Optional string. The obs column holding the clustering. Defaults to
  the one the PAGA run used, and errors if a different one is given:
  nodes from one clustering with edges from another produce a plausible
  looking graph that is wrong.

- threshold:

  Numeric. Edges below this connectivity are dropped. Defaults to
  `0.01`.

- tree_only:

  Boolean. Use the maximum spanning forest rather than the full
  abstracted graph. Defaults to `FALSE`.

- centroid:

  String. One of `c("median", "mean")`. How a cluster's position is
  summarised. Median by default, since embeddings throw stragglers that
  drag a mean off its cluster.

- node_stat_col:

  Optional string. A numeric obs column to summarise per cluster with
  the same statistic, e.g. `"palantir_pseudotime"`, giving a `stat`
  column to colour the nodes by.

- ...:

  Additional arguments forwarded to
  [`extract_embedding_data()`](https://gregorlueg.github.io/bixverse/reference/extract_embedding_data.md)
  and onward to
  [`get_embedding()`](https://gregorlueg.github.io/bixverse/reference/get_embedding.md)
  (e.g. `modality`).

## Value

A list with the embedding stored as an `embedding` attribute and

- nodes - data.table with `cluster` (a factor in graph order), `dim_1`,
  `dim_2`, `n_cells` and, when `node_stat_col` is given, `stat`.
  Clusters holding no cells are dropped, as they have no position.

- edges - data.table with `from`, `to`, `weight` and the coordinates
  `x`, `y`, `xend`, `yend` of both ends, ready for a segment layer. Each
  edge appears once, despite the graph being stored symmetrically.

## References

Wolf, et al., Genome Biol., 2019.
