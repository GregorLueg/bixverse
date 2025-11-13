# Get correlation-based graph

Helper function to get a correlation-based igraph from the class

## Usage

``` r
get_cor_graph(object, epsilon, .verbose)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md).

- epsilon:

  Float. The epsilon parameter for the RBF function, in this case the
  bump function.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the following elements:

- graph - The igraph

- params - A list that contains the parameters of the graph generation
  and general graph information (node, edge numbers).
