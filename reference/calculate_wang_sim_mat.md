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

## Examples

``` r
# Wang similarity matrix with relationship-specific weights
onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f"),
  type = c("part_of", "part_of", "part_of", "is_a", "is_a")
)
weights <- c("part_of" = 0.8, "is_a" = 0.6)
round(calculate_wang_sim_mat(onto, weights = weights), 3)
#>       f     d     b     a     c     e
#> f 1.000 0.470 0.625 0.400 0.796 0.428
#> d 0.470 1.000 0.764 0.477 0.590 0.558
#> b 0.625 0.764 1.000 0.643 0.764 0.742
#> a 0.400 0.477 0.643 1.000 0.477 0.481
#> c 0.796 0.590 0.764 0.477 1.000 0.558
#> e 0.428 0.558 0.742 0.481 0.558 1.000
```
