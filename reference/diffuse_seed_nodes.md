# Diffuse seed genes over a network

This function takes a diffusion vector and leverages personalised
page-rank diffusion to identify influential nodes. These can be used
subsequently for community detection or check AUROC values given a set
of genes.

## Usage

``` r
diffuse_seed_nodes(
  object,
  diffusion_vector,
  summarisation = c("max", "mean", "harmonic_sum")
)
```

## Arguments

- object:

  `network_diffusions` object. The underlying class
  [`network_diffusions()`](network_diffusions.md).

- diffusion_vector:

  Named nuermic. A named vector with values to use for the reset
  parameter in the personalised page-rank diffusion. Names should
  represent node names of the graph.

- summarisation:

  String. If there are duplicated names in the `diffusion_vector` how to
  summarise the scores.

## Value

The class with added diffusion score based on a single set of seed
genes. Additionally, the seed genes are stored in the class.
