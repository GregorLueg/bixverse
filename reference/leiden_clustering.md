# Leiden or Louvain community detection on an igraph

Modularity-based Leiden (Traag et al.) or Louvain (Blondel et al.)
community detection, run in Rust. Mirrors
`igraph::cluster_leiden(objective_function = "modularity")`, but the
result depends on `seed` only (not on the global RNG or the number of
threads). Communities are numbered by size, `1` being the largest.

## Usage

``` r
leiden_clustering(
  graph,
  resolution = 1,
  algorithm = c("leiden", "louvain"),
  beta = 0.01,
  n_iterations = 2L,
  seed = 42L
)
```

## Arguments

- graph:

  igraph. Undirected. If it has an edge attribute `weight`, it is used,
  otherwise all edges get weight 1.

- resolution:

  Numeric. Modularity resolution parameter. Higher values give more
  communities.

- algorithm:

  String. One of `c("leiden", "louvain")`.

- beta:

  Numeric. Randomness of the Leiden refinement merges. Ignored for
  Louvain.

- n_iterations:

  Integer. Number of full iterations, each starting from the previous
  partition.

- seed:

  Integer. Random seed.

## Value

Integer vector with the community of each node (1-based, `1` is the
largest community).

## References

Traag, Waltman & van Eck, Sci Rep, 2019; Blondel, et al., J Stat Mech,
2008

## Examples

``` r
g <- igraph::make_graph("Zachary")
table(leiden_clustering(g))
#> 
#>  1  2  3  4 
#> 12 11  6  5 
```
