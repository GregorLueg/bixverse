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
  n_iterations = NULL,
  min_gain = 0,
  local_moving = c("queue", "coloured"),
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

  Optional integer. Number of full iterations, each starting from the
  previous partition. A negative value iterates until the partition is
  stable. If `NULL`, `5L` for Leiden and `2L` for Louvain, which gave
  the best modularity across synthetic graphs.

- min_gain:

  Numeric. With a negative `n_iterations`, also stop once an iteration
  gains no more than this in modularity.

- local_moving:

  String. One of `c("queue", "coloured")`. `"queue"` is Traag's
  sequential queue. `"coloured"` moves the nodes of one colour class in
  parallel, which only kicks in on larger levels of the graph.

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
