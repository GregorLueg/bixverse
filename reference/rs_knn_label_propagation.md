# kNN label propagation

The function is a helper function to do kNN label propagation. This can
be useful for semi-supervised tasks. It implements the label spreading
method.

## Usage

``` r
rs_knn_label_propagation(
  edge_list,
  one_hot_encoding,
  label_mask,
  alpha,
  iterations,
  tolerance
)
```

## Arguments

- edge_list:

  Integer vector. In form of node_1, node_2, node_3, ... which indicates
  alternating pairs (node_1, node_2), etc in terms of edges

- one_hot_encoding:

  Integer matrix. Each row represents a sample, the columns the one-hot
  encodings. Everything 0 denotes the unlabelled data.

- label_mask:

  Boolean vector. Which of the samples do not have a label. Needs to be
  same length as `nrow(one_hot_encoding)`.

- alpha:

  Numeric. Parameter that controls the spreading. Usually between 0.9 to
  0.95. Larger values drive further labelling, smaller values are more
  conversative.

- iterations:

  For how many (max) iterations to run the algorithm.

- tolerance:

  If the value below this is reached, an early stop is initialised

## Value

The matrix with the probabilities of being of a certain class
