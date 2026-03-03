# kNN label propagation

The function is a helper function to do kNN label propagation. This can
be useful for semi-supervised tasks. It implements the label spreading
method.

## Usage

``` r
rs_knn_label_propagation(
  from,
  to,
  one_hot_encoding,
  label_mask,
  weights,
  label_prop_params
)
```

## Arguments

- from:

  Integer vector. Source node indices for each edge.

- to:

  Integer vector. Target node indices for each edge. Must be the same
  length as `from`.

- one_hot_encoding:

  Integer matrix. Each row represents a sample, the columns the one-hot
  encodings. Everything 0 denotes the unlabelled data.

- label_mask:

  Boolean vector. Which of the samples do not have a label. Needs to be
  same length as `nrow(one_hot_encoding)`.

- weights:

  Optional numeric vector. Edge weights for each pair in `from`/`to`.
  Must have the same length as `from`. If NULL, all edges are treated as
  unweighted.

- label_prop_params:

  List. Named list of parameters with the following optional fields
  (defaults in parentheses):

  - `alpha` numeric, spreading strength (0.9)

  - `iter` integer, max iterations (100)

  - `tolerance` numeric, convergence threshold (1e-6)

  - `symmetrise` logical, symmetrise the graph (FALSE)

  - `symmetry_strategy` character, one of "average", "min", "max"
    ("average")

  - `max_hops` integer, restrict spreading radius (unrestricted)

## Value

The matrix with the probabilities of being of a certain class.
