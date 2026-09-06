# Get the diffusion vector

Returns the diffusion vector if you ran
[`tied_diffusion()`](https://gregorlueg.github.io/bixverse/reference/tied_diffusion.md)
or
[`diffuse_seed_nodes()`](https://gregorlueg.github.io/bixverse/reference/diffuse_seed_nodes.md).

## Usage

``` r
get_diffusion_vector(object)
```

## Arguments

- object:

  The underlying class
  [`NetworkDiffusions()`](https://gregorlueg.github.io/bixverse/reference/NetworkDiffusions.md).

## Value

The diffusion vector if found. If you did not run either diffusion
functions, it will return `NULL` and a warning.

## Examples

``` r
# diffusion scores after a single seed node diffusion
set.seed(42)
g <- igraph::sample_pa(15, directed = FALSE)
edges <- data.table::setDT(igraph::as_data_frame(g))[, `:=`(
  from = sprintf("node_%i", from),
  to = sprintf("node_%i", to)
)]
object <- NetworkDiffusions(edges, weighted = FALSE, directed = FALSE)
object <- diffuse_seed_nodes(object, c(node_1 = 1, node_3 = 1), "max")
head(get_diffusion_vector(object))
#>     node_1     node_2     node_3     node_4     node_6     node_7 
#> 0.12792857 0.24907561 0.15054927 0.10645037 0.04364123 0.03541402 
```
