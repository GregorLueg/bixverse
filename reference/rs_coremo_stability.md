# Helper function to assess CoReMo cluster stability

This function is a helper for the leave-on-out stability assessment of
CoReMo clusters. The function will generate the distance vectors based
on leaving out the samples defined in indices one by one.

## Usage

``` r
rs_coremo_stability(data, indices, epsilon, rbf_type, spearman)
```

## Arguments

- data:

  Numeric matrix. The original processed matrix.

- indices:

  Integer vector. The sample indices to remove to re-calculate the
  distances.

- epsilon:

  Float. Epsilon parameter for the RBF.

- rbf_type:

  String. Needs to be from `c("gaussian", "bump", "inverse_quadratic")`.

- spearman:

  Boolean. Shall Spearman correlation be used.

## Value

A list with `length(indices)` elements, each containing the distance
minus the given sample.
