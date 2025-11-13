# Calculate massively parallelised personalised page rank scores

Helper function to calculate in parallel on the same (unweighted)
network the personalised page rank as fast as possible. Can be used for
permutations type approaches.

## Usage

``` r
rs_page_rank_parallel(
  node_names,
  from,
  to,
  weights,
  diffusion_scores,
  undirected
)
```

## Arguments

- node_names:

  String vector. Name of the graph nodes.

- from:

  String vector. The names of the `from` edges from the edge list.

- to:

  String vector. The names of the `to` edges from the edge list.

- weights:

  Optional weight vector. If NULL, defaults to 1.0 as weight for all
  edges.

- diffusion_scores:

  List. The personalised vectors for the page rank reset values. Each
  element must sum to 1 and be of same length of `node_names`!

- undirected:

  Boolean. Is this an undirected graph.

## Value

A matrix of the scores with each row representing an element in the
`diffusion_scores` list (in order), and each column representing the
value of the personalised page rank diffusion for this node.
