# Helper to get the Gene Ontology levels

Gets the ontology depth based on the edge dt

## Usage

``` r
get_go_levels(edge_dt)
```

## Arguments

- edge_dt:

  data.table. The gene ontology edge data (i.e., term connections
  between the different terms)

## Value

A data.table with the identifier and depth.
