# Calculate the Wang similarities between terms

This function calculates the Wang similarity for a set of temrs, based
on the DAG for a given ontology.

## Usage

``` r
calculate_wang_sim(terms, parent_child_dt, weights, add_self = FALSE)
```

## Arguments

- terms:

  String vector. The terms for which to calculate the Wang similarity.

- parent_child_dt:

  data.table. The data.table with column parent and child. You also need
  to have a type column for the Wang similarity to provide the weights
  for the relationships.

- weights:

  Named numeric. The relationship of type to weight for this specific
  edge. For example `c("part_of" = 0.8, "is_a" = 0.6)`.

- add_self:

  Boolean. Shall self-similarities be added. Defaults to `FALSE`.

## Value

A data.table with the calculated similarities.
