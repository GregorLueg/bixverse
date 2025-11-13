# Apply a Radial Basis Function

Applies a radial basis function (RBF) to a given distance vector. Has at
the option to apply a Gaussian, Bump or Inverse Quadratic RBF.

## Usage

``` r
rs_rbf_function(x, epsilon, rbf_type)
```

## Arguments

- x:

  Numeric vector. The distances you wish to apply the Gaussian kernel
  onto.

- epsilon:

  Float. Epsilon parameter for the RBF.

- rbf_type:

  String. Needs to be from `c("gaussian", "bump", "inverse_quadratic")`.

## Value

The affinities after the Kernel was applied.
