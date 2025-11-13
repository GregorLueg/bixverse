# kNN-based graph label propagation

In case of a kNN graph with a subset of nodes having labels, this
function can be used to propagate the labels through the graph. This
function uses the label spreading version of the algorithm.

## Usage

``` r
knn_graph_label_propagation(
  edge_list,
  labels,
  iterations = 100L,
  alpha = 0.9,
  tolerance = 1e-06
)
```

## Arguments

- edge_list:

  Integer vector. In form of node_1, node_2, node_3, ... which indicates
  alternating pairs (node_1, node_2), etc in terms of edges.

- labels:

  String. The labels with `NA` indicating unlabelled.

- iterations:

  Integer. Maximum iterations for the algorithm.

- alpha:

  Numeric. Parameter that controls the spreading. Usually between 0.9 to
  0.95. Larger values drive further labelling, smaller values are more
  conversative.

- tolerance:

  Tolerance parameter for early stopping.

## Value

A list with the following elements

- assignment_probs - Matrix with the assignment probabilities.

- final_labels - Final labels in the graph
