# Helper function to assess CoReMo cluster stability

**\[experimental\]** This function is a helper for the leave-one-out
stability assessment of CoReMo clusters. The function will generate the
distance vectors based on leaving out the samples defined in indices one
by one. Distances are `1 - rbf(1 - |cor|)` over the feature
correlations.

## Usage

``` r
rs_coremo_stability(data, indices, epsilon, rbf_type, spearman)
```

## Arguments

- data:

  Numeric matrix. The original processed matrix, samples x features.

- indices:

  Integer vector. The 1-based sample (row) indices to remove, one at a
  time, to re-calculate the distances.

- epsilon:

  Float. Epsilon parameter for the RBF.

- rbf_type:

  String. Needs to be from `c("gaussian", "bump", "inverse_quadratic")`.

- spearman:

  Boolean. Shall Spearman correlation be used.

## Value

A list with `length(indices)` elements, each containing the flattened
upper-triangle feature distances with that sample removed.
