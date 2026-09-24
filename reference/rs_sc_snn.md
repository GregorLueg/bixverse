# Generates the sNN graph for igraph

**\[experimental\]** This function takes a kNN matrix and generates the
inputs for an SNN graph based on it.

## Usage

``` r
rs_sc_snn(knn_mat, snn_method, limited_graph, pruning, verbose)
```

## Arguments

- knn_mat:

  Integer matrix. Rows represent cells and the columns represent the
  neighbours (0-indexed!).

- snn_method:

  String. Which method to use to calculate the similarity. Choice of
  `c("jaccard", "rank")`; anything else errors.

- limited_graph:

  Boolean. Shall the sNNs only be calculated between direct neighbours
  in the graph, or between all possible combinations.

- pruning:

  Float. Below which similarity value to prune the weight to 0.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items:

- edges - Integer vector with the flattened sNN edge pairs (1-indexed!),
  ready for
  [`igraph::add_edges()`](https://r.igraph.org/reference/add_edges.html).

- weights - sNN weights of the pairs above.
