# Constrained personalised page rank over a list

This function implements a personalised constrained page-rank over a
list of personalisation vectors for the same network. You can define
`sink_nodes` or `sink_edges`. In the case of the former, the vist of a
sink_node (via a type attribute in the igraph) will automatically cause
the surfer to reset. In the case of the latter, the traversal of a
sink_edge is allowed, however, subsequently, the surfer will be reseted.
This function can be useful for deriving PPR profiles under constraints
in heterogenous graphs and has been inspired by Ruiz, et al.

## Usage

``` r
constrained_page_rank_ls(
  graph,
  personalisation_list,
  sink_nodes = NULL,
  sink_edges = NULL
)
```

## Arguments

- graph:

  igraph. This one needs to be directed and weighted and have the node
  attribute `type` defining the node type and the edge attribute `type`
  defining the edge type.

- personalisation_list:

  A list of numerical vectors to use for the personalisation.

- sink_nodes:

  Optional String vector. The node types that should force a reset.

- sink_edges:

  Optional String vector. The edge types after which there should be a
  forced reset.

## Value

A list with the constrained personalised page rank values.

## References

Ruiz, et al., Nat Commun, 2021

## Examples

``` r
# two personalisation vectors over the same signalling graph
nodes <- data.frame(
  name = c("rec_a", "kin_b", "tf_c", "gene_d"),
  type = c("receptor", "kinase", "tf", "target_gene")
)
edges <- data.frame(
  from = c("rec_a", "kin_b", "tf_c"),
  to = c("kin_b", "tf_c", "gene_d"),
  weight = rep(1, 3),
  type = c("activation", "phosphorylation", "tf_activation")
)
g <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
res <- constrained_page_rank_ls(
  g,
  personalisation_list = list(a = c(1, 0, 0, 0), b = c(0, 1, 0, 0))
)
res$a
#>     rec_a     kin_b      tf_c    gene_d 
#> 0.3138116 0.2667399 0.2267289 0.1927196 
```
