# Return ancestry terms from an ontology

This function will return all ancestors and descendants based on a
provided data.table with parent-child terms

## Usage

``` r
get_ontology_ancestry(parent_child_dt)
```

## Arguments

- parent_child_dt:

  data.table. The data.table with column parent and child.

## Value

A list with

- ancestors A list with all ancestor terms.

- descendants A list with all descendant terms.

## Examples

``` r
# ancestors and descendants of every term in a toy ontology
onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f")
)
ancestry <- get_ontology_ancestry(onto)
ancestry$ancestors[["f"]]
#> [1] "f" "c" "b" "a"
```
