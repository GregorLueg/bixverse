# Apply a Radial Basis Function

**\[experimental\]** Applies a radial basis function (RBF) to a given
distance vector. Has the option to apply a Gaussian, Bump or Inverse
Quadratic RBF.

## Usage

``` r
rs_rbf_function(x, epsilon, rbf_type)
```

## Arguments

- x:

  Numeric vector. The distances you wish to apply the RBF onto.

- epsilon:

  Float. Epsilon parameter for the RBF.

- rbf_type:

  String. Needs to be from `c("gaussian", "bump", "inverse_quadratic")`.
  Other values raise an error.

## Value

The affinities after the Kernel was applied.
