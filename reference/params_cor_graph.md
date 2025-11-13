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

  Float. Defines the epsilon parameter for the radial basis function.
  Defaults to 2, but should be ideally optimised.

- min_cor:

  Float. Minimum absolute correlation that needs to be observed in
  either data set. Only relevant for differential correlation-based
  graphs.

- fdr_threshold:

  Float. Maximum FDR for the differential correlation p-value.

- verbose:

  Boolean. Controls verbosity of the graph generation function.

## Value

List with parameters for usage in subsequent function.
