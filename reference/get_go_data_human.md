# Wrapper function to load and process the gene ontology data.

This function loads in gene ontology data stored in the package and
processes it into the format for
[`gene_ontology_data()`](gene_ontology_data.md). Wraps
[`load_go_human_data()`](load_go_human_data.md) and
[`process_go_data()`](process_go_data.md) into one.

## Usage

``` r
get_go_data_human(filter_relationships = TRUE, .verbose = TRUE)
```

## Arguments

- filter_relationships:

  Boolean. Shall the ontology be filtered to only `"is_a"` and
  `"part_of"` relationships. Defaults to TRUE.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A data.table
