# Find RBH communities

This function will identify communities in the reciprocal best hit (RBH)
graph. It will iterate through resolutions and add the results to the
class. Additionally, a column will be added that signifies the
resolution with the best modularity.

## Usage

``` r
find_rbh_communities(
  object,
  resolution_params = params_graph_resolution(),
  max_workers = NULL,
  parallel = TRUE,
  random_seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  The underlying class, see
  [`RbhGraph()`](https://gregorlueg.github.io/bixverse/reference/RbhGraph.md).

- resolution_params:

  List. Parameters for the resolution search, see
  [`params_graph_resolution()`](https://gregorlueg.github.io/bixverse/reference/params_graph_resolution.md).
  Contains:

  - min_res - Float. Minimum resolution to test.

  - max_res - Float. Maximum resolution to test.

  - number_res - Integer. Number of resolutions to test between the
    `max_res` and `min_res.`

- max_workers:

  Integer. Number of maximum cores to use. Defaults to half of the
  identified cores (to a maximum of 8).

- parallel:

  Boolean. Shall the resolution search be in parallel.

- random_seed:

  Integer. Random seed for reproducibility.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added community detection results.

## Examples

``` r
# Leiden communities across a resolution sweep of the RBH graph
set.seed(123)
modules <- data.table::data.table(
  origin = rep(c("set_a", "set_b"), each = 20),
  module = rep(c("m1", "m2", "m3", "m4"), each = 10),
  gene = unlist(replicate(4, sample(letters, 10), simplify = FALSE))
)
object <- RbhGraph(
  modules,
  rbh_type = "set",
  dataset_col = "origin",
  module_col = "module",
  value_col = "gene"
)
object <- generate_rbh_graph(object, minimum_similarity = 0)
object <- find_rbh_communities(object, parallel = FALSE, .verbose = FALSE)
head(get_results(object))
#>    node_name origin_id module_id membership resolution modularity
#>       <char>    <char>    <char>      <num>      <num>      <num>
#> 1:  set_a_m1     set_a        m1          1  0.1000000  0.4897959
#> 2:  set_a_m2     set_a        m2          2  0.1000000  0.4897959
#> 3:  set_b_m3     set_b        m3          1  0.1000000  0.4897959
#> 4:  set_b_m4     set_b        m4          2  0.1000000  0.4897959
#> 5:  set_a_m1     set_a        m1          1  0.1389495  0.4897959
#> 6:  set_a_m2     set_a        m2          2  0.1389495  0.4897959
#>    best_modularity
#>             <lgcl>
#> 1:            TRUE
#> 2:            TRUE
#> 3:            TRUE
#> 4:            TRUE
#> 5:            TRUE
#> 6:            TRUE
```
