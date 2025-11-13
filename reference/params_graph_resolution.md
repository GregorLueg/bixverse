# Wrapper function to generate resolution parameters for Leiden or Louvain clustering.

Wrapper function to generate resolution parameters for Leiden or Louvain
clustering.

## Usage

``` r
params_graph_resolution(min_res = 0.1, max_res = 10, number_res = 15L)
```

## Arguments

- min_res:

  Float. Minimum resolution to test.

- max_res:

  Float. Maximum resolution to test.

- number_res:

  Integer. Number of resolutions to test between the `max_res` and
  `min_res.`

## Value

List with parameters for usage in subsequent function.
