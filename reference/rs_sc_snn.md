# Generates the sNN graph for igraph

This function takes a kNN matrix and generates the inputs for an SNN
graph based on it.

## Usage

``` r
rs_sc_snn(knn_mat, snn_method, limited_graph, pruning, verbose)
```

## Arguments

- knn_mat:

  Integer matrix. Rows represent cells and the columns represent the
  neighbours.

- snn_method:

  String. Which method to use to calculate the similarity. Choice of
  `c("jaccard", "rank")`.

- limited_graph:

  Boolean. Shall the sNNs only be calculated between direct neighbours
  in the graph, or between all possible combinations.

- pruning:

  Float. Below which value for the Jaccard similarity to prune the
  weight to 0.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the following items:

- edges - sNN edges as edge pairs.

- weights - sNN weights of the pairs above.
