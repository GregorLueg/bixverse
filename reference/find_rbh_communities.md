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

  The underlying class, see [`rbh_graph()`](rbh_graph.md).

- resolution_params:

  List. Parameters for the resolution search, see
  [`params_graph_resolution()`](params_graph_resolution.md). Contains:

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
