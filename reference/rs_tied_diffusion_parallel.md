# Calculate massively parallelised tied diffusion scores

Helper function to calculate in parallel on the same (unweighted)
network the tied diffusions as fast as possible. Can be used for
permutation.

## Usage

``` r
rs_tied_diffusion_parallel(
  node_names,
  from,
  to,
  weights,
  diffusion_scores_1,
  diffusion_scores_2,
  summarisation_fun,
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

- diffusion_scores_1:

  List. The first set of personalised vectors for the page rank reset
  values. Each element must sum to 1 and be of same length of
  `node_names`!

- diffusion_scores_2:

  List. The second set of personalised vectors for the page rank reset
  values. Each element must sum to 1 and be of same length of
  `node_names`!

- summarisation_fun:

  String. One of `c("min", "max", "avg")`. Which type of summarisation
  function to use to calculate the tied diffusion.

- undirected:

  Boolean. Is this an undirected graph.

## Value

A matrix of the scores with each row representing a tied diffusion of of
`diffusion_scores_1` and `diffusion_scores_2` lists (in order), and each
column representing the value of the tied diffusion for this node.
