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

  `NetworkDiffusions` object. The underlying class
  [`NetworkDiffusions()`](https://gregorlueg.github.io/bixverse/reference/NetworkDiffusions.md).

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

## Examples

``` r
# personalised page-rank diffusion from three seed nodes
set.seed(42)
g <- igraph::sample_pa(15, directed = FALSE)
edges <- data.table::setDT(igraph::as_data_frame(g))[, `:=`(
  from = sprintf("node_%i", from),
  to = sprintf("node_%i", to)
)]
object <- NetworkDiffusions(edges, weighted = FALSE, directed = FALSE)
object <- diffuse_seed_nodes(
  object,
  c(node_1 = 1, node_3 = 1, node_10 = 1),
  summarisation = "max"
)
head(get_diffusion_vector(object))
#>     node_1     node_2     node_3     node_4     node_6     node_7 
#> 0.10750018 0.27058907 0.12650867 0.08945174 0.03667234 0.02975890 
```
