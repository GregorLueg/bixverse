# Leiden or Louvain community detection

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
  local_moving = NULL,
  node_order = NULL,
  seed = 42L
)
```

## Arguments

- graph:

  `SingleCellConnectivity` or undirected igraph. For an igraph, the edge
  attribute `weight` is used if present, otherwise all edges get
  weight 1. A `SingleCellConnectivity` is clustered in Rust directly; an
  igraph is flattened into an edge list first. Both give identical
  labels for the same graph.

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
  stable. If `NULL`, `2L`.

- min_gain:

  Numeric. With a negative `n_iterations`, also stop once an iteration
  gains no more than this in modularity.

- local_moving:

  Optional string. One of `c("queue", "coloured")`. `"queue"` is Traag's
  sequential queue. `"coloured"` moves the nodes of one colour class in
  parallel, which only kicks in on larger levels of the graph. `NULL` is
  `"queue"` for Leiden and `"coloured"` for Louvain.

- node_order:

  Optional string. Visit order of the queue, `"bfs"` (sweep in BFS order
  from a random start) or `"shuffle"` (as igraph). `NULL` is `"bfs"` for
  Leiden and `"shuffle"` for Louvain.

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
