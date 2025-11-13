# Helper to identify the right epsilon parameter

This function will take a distance vector from the upper triangle of a
symmetric distance matrix and apply the desired RBF with the supplied
epsilon from epsilon vec. Subsequently, the column sums will be measured
to identify the total similarity of each feature with other features.
This data can be used to see if the data follows scale-free topology for
example to identify the right epsilon parameter with the given RBF.

## Usage

``` r
rs_rbf_iterate_epsilons(dist, epsilon_vec, original_dim, shift, rbf_type)
```

## Arguments

- dist:

  Numeric vector. The distances you wish to apply the RBF function to.

- epsilon_vec:

  Numeric vector. The epsilons you wish to use/test.

- original_dim:

  Integer. The original dimensions of the symmetric distance matrix.

- shift:

  Integer. Was the matrix shifted up (0 = diagonal included; 1 diagonal
  not incldued).

- rbf_type:

  String. Option of `c('gaussian', 'bump')` for the currently
  implemented RBF function.

## Value

A matrix with rows being the epsilons tested, and columns representing
the summed affinity to other features.
