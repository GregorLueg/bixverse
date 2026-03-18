# Ontology class (deprecated)

**\[deprecated\]**

This constructor has been renamed to [`OntologySim()`](OntologySim.md).

## Usage

``` r
ontology(parent_child_dt, .verbose = TRUE)
```

## Arguments

- parent_child_dt:

  A data.table containing the ontological information as parent-child
  relationships. Must contain the `c("parent", "child")` columns.

- .verbose:

  Boolean. Controls the verbosity of the class.

## Value

Returns a [`OntologySim()`](OntologySim.md) object.
