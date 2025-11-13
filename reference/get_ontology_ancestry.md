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
