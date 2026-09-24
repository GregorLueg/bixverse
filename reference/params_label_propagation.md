# Wrapper function to generate label propagation parameters

Wrapper function to generate label propagation parameters

## Usage

``` r
params_label_propagation(
  alpha = 0.9,
  iter = 100L,
  tolerance = 1e-06,
  symmetrise = TRUE,
  symmetry_strategy = c("average", "avg", "min", "max"),
  max_hops = NULL
)
```

## Arguments

- alpha:

  Numeric. Controls the spreading strength. Higher values anchor
  labelled nodes more strongly to their original label. Defaults to
  `0.9`.

- iter:

  Integer. Maximum number of iterations to run. Defaults to `100L`.

- tolerance:

  Numeric. Convergence threshold. Stops early if the maximum change
  across all nodes falls below this value. Defaults to `1e-06`.

- symmetrise:

  Boolean. Whether to symmetrise the graph. Defaults to `TRUE`.

- symmetry_strategy:

  String. Strategy to resolve weight conflicts when symmetrising. Only
  relevant when `symmetrise = TRUE` and edge weights are provided. One
  of `c("average", "avg", "min", "max")`. Defaults to `"average"`.

- max_hops:

  Integer or `NULL`. If provided, restricts label spreading to nodes
  within this many hops of any labelled node. Nodes beyond this limit
  remain all-zeroes. Defaults to `NULL`.

## Value

A named list with the following elements:

- alpha - Numeric. Controls the spreading strength. Higher values anchor
  labelled nodes more strongly to their original label. Defaults to
  `0.9`.

- iter - Integer. Maximum number of iterations to run. Defaults to
  `100L`.

- tolerance - Numeric. Convergence threshold. Stops early if the maximum
  change across all nodes falls below this value. Defaults to `1e-06`.

- symmetrise - Boolean. Whether to symmetrise the graph. Defaults to
  `TRUE`.

- symmetry_strategy - String. Strategy to resolve weight conflicts when
  symmetrising. Only relevant when `symmetrise = TRUE` and edge weights
  are provided. One of `c("average", "avg", "min", "max")`. Defaults to
  `"average"`.

- max_hops - Integer or `NULL`. If provided, restricts label spreading
  to nodes within this many hops of any labelled node. Nodes beyond this
  limit remain all-zeroes. Defaults to `NULL`.
