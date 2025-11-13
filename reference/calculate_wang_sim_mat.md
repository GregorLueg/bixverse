# Calculate the Wang similarity matrix

This function calculates the Wang similarity, based on the DAG for a
given ontology. This function will return the full similarity matrix.

## Usage

``` r
calculate_wang_sim_mat(parent_child_dt, weights)
```

## Arguments

- parent_child_dt:

  data.table. The data.table with column parent and child. You also need
  to have a type column for the Wang similarity to provide the weights
  for the relationships.

- weights:

  Named numeric. The relationship of type to weight for this specific
  edge. For example `c("part_of" = 0.8, "is_a" = 0.6)`.

## Value

The symmetric Wang similarity matrix.
