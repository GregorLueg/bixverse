# Calculate the AUROC for a diffusion score

This functions can take a given `NetworkDiffusions` object and
calculates an AUC and generates a Z-score based on random permutation of
`random_aucs` for test for statistical significance if desired.

## Usage

``` r
calculate_diffusion_auc(
  object,
  hit_nodes,
  auc_iters = 10000L,
  random_aucs = 1000L,
  permutation_test = FALSE,
  seed = 42L
)
```

## Arguments

- object:

  `NetworkDiffusions` object. The underlying class
  [`NetworkDiffusions()`](https://gregorlueg.github.io/bixverse/reference/NetworkDiffusions.md).

- hit_nodes:

  String vector. Which nodes in the graph are considered a 'hit'.

- auc_iters:

  Integer. How many iterations to run to approximate the AUROC.

- random_aucs:

  Integer. How many random AUROCs to calculate to estimate the Z-score.
  Only of relevance if permutation test is set to `TRUE`.

- permutation_test:

  Boolean. Shall a permutation based Z-score be calculated.

- seed:

  Integer. Random seed.

## Value

List with AUC and Z-score as the two named elements if permutations test
set to TRUE; otherwise just the AUC.

## Examples

``` r
# AUROC of the diffusion score against two known hit nodes
set.seed(42)
g <- igraph::sample_pa(15, directed = FALSE)
edges <- data.table::setDT(igraph::as_data_frame(g))[, `:=`(
  from = sprintf("node_%i", from),
  to = sprintf("node_%i", to)
)]
object <- NetworkDiffusions(edges, weighted = FALSE, directed = FALSE)
object <- diffuse_seed_nodes(object, c(node_1 = 1, node_3 = 1), "max")
calculate_diffusion_auc(object, hit_nodes = c("node_2", "node_4"))
#> [1] 0.9239
```
