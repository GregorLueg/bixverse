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

## Examples

``` r
# Wang similarity for a subset of terms only
onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f"),
  type = c("part_of", "part_of", "part_of", "is_a", "is_a")
)
weights <- c("part_of" = 0.8, "is_a" = 0.6)
calculate_wang_sim(c("c", "d", "f"), onto, weights = weights)
#>     term1  term2      sims
#>    <char> <char>     <num>
#> 1:      c      d 0.5901639
#> 2:      c      f 0.7960848
#> 3:      d      f 0.4698206
```
