# Wrapper function to generate resolution parameters for Leiden or Louvain clustering.

Wrapper function to generate resolution parameters for Leiden or Louvain
clustering.

## Usage

``` r
params_graph_resolution(min_res = 0.1, max_res = 10, number_res = 15L)
```

## Arguments

- min_res:

  Numeric. Minimum resolution to test. Defaults to `0.1`.

- max_res:

  Numeric. Maximum resolution to test. Defaults to `10.0`.

- number_res:

  Integer. Number of resolutions to test between the `max_res` and
  `min_res.` Defaults to `15L`.

## Value

A named list with the following elements:

- min_res - Numeric. Minimum resolution to test. Defaults to `0.1`.

- max_res - Numeric. Maximum resolution to test. Defaults to `10.0`.

- number_res - Integer. Number of resolutions to test between the
  `max_res` and `min_res.` Defaults to `15L`.
