# Helper function to plot distance to affinity relationship

This function plots the distance to affinity relationship after applying
a given radial basis function (RBF) with a given epsilon value.

## Usage

``` r
plot_rbf_impact(rbf_type, epsilon)
```

## Arguments

- rbf_type:

  String. The RBF type. One of
  `c("gaussian", "bump", "inverse_quadratic")`

- epsilon:

  Numerical. The epsilon parameter.

## Value

A plot depicting the distance to affinity relationship after applying
the RBF function.
