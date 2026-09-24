# Wrapper function for graph generation

Wrapper function for graph generation

## Usage

``` r
params_cor_graph(
  epsilon = 2,
  min_cor = 0.2,
  fdr_threshold = 0.05,
  verbose = TRUE
)
```

## Arguments

- epsilon:

  Numeric. Defines the epsilon parameter for the radial basis function.
  Defaults to `2.0`.

- min_cor:

  Numeric. Minimum absolute correlation that needs to be observed in
  either data set. Only relevant for differential correlation-based
  graphs. Defaults to `0.2`.

- fdr_threshold:

  Numeric. Maximum FDR for the differential correlation p-value.
  Defaults to `0.05`.

- verbose:

  Boolean. Controls verbosity of the graph generation function. Defaults
  to `TRUE`.

## Value

A named list with the following elements:

- epsilon - Numeric. Defines the epsilon parameter for the radial basis
  function. Defaults to `2.0`.

- min_cor - Numeric. Minimum absolute correlation that needs to be
  observed in either data set. Only relevant for differential
  correlation-based graphs. Defaults to `0.2`.

- fdr_threshold - Numeric. Maximum FDR for the differential correlation
  p-value. Defaults to `0.05`.

- verbose - Boolean. Controls verbosity of the graph generation
  function. Defaults to `TRUE`.
