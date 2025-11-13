# Diffuse seed genes in a tied manner over a network

This function takes two sets of diffusion vector and leverages tied
diffusion to identify an intersection of influential nodes. If the
network is undirected, the method will run two personalised page rank
diffusions based on the diffusion vectors and generate the score
aggregation

## Usage

``` r
tied_diffusion(
  object,
  diffusion_vector_1,
  diffusion_vector_2,
  summarisation = c("max", "mean", "harmonic_sum"),
  score_aggregation = c("min", "max", "mean"),
  .verbose = FALSE
)
```

## Arguments

- object:

  `network_diffusions` object. The underlying class
  [`network_diffusions()`](network_diffusions.md).

- diffusion_vector_1:

  Named numeric. The first named vector with values to use for the reset
  parameter in the personalised page-rank diffusion. Names should
  represent node names of the graph.

- diffusion_vector_2:

  Named numeric. The second named vector with values to use for the
  reset parameter in the personalised page-rank diffusion. Names should
  represent node names of the graph.

- summarisation:

  String. If there are duplicated names in the `diffusion_vector` how to
  summarise these.

- score_aggregation:

  String. How to summarise the tied scores.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added diffusion score based on a two sets of seed genes.
Additionally, the seed genes are stored in the class.
