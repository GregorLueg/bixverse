# Apply a Radial Basis Function (to a matrix)

Applies a radial basis function (RBF) to a given distance matrix. Has at
the option to apply a Gaussian, Bump or Inverse Quadratic RBF.

## Usage

``` r
rs_rbf_function_mat(x, epsilon, rbf_type)
```

## Arguments

- x:

  Numeric Matrix. The distances you wish to apply the Gaussian kernel
  onto.

- epsilon:

  Float. Epsilon parameter for the RBF.

- rbf_type:

  String. Needs to be from `c("gaussian", "bump", "inverse_quadratic")`.

## Value

The affinities after the Kernel was applied.
