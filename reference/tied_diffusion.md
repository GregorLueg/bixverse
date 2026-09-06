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

  `NetworkDiffusions` object. The underlying class
  [`NetworkDiffusions()`](https://gregorlueg.github.io/bixverse/reference/NetworkDiffusions.md).

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

## Examples

``` r
# tied diffusion between two sets of seed nodes
set.seed(42)
g <- igraph::sample_pa(15, directed = FALSE)
edges <- data.table::setDT(igraph::as_data_frame(g))[, `:=`(
  from = sprintf("node_%i", from),
  to = sprintf("node_%i", to)
)]
object <- NetworkDiffusions(edges, weighted = FALSE, directed = FALSE)
object <- tied_diffusion(
  object,
  diffusion_vector_1 = c(node_1 = 1, node_3 = 1),
  diffusion_vector_2 = c(node_2 = 1, node_6 = 1),
  summarisation = "max",
  score_aggregation = "min"
)
head(get_diffusion_vector(object))
#>     node_1     node_2     node_3     node_4     node_6     node_7 
#> 0.04299109 0.20231102 0.06764978 0.10645037 0.04364123 0.03541402 
```
