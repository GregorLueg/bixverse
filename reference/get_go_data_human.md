# Wrapper function to load and process the gene ontology data.

This function loads in gene ontology data stored in the package and
processes it into the format for
[`GeneOntologyElim()`](https://gregorlueg.github.io/bixverse/reference/GeneOntologyElim.md).
Wraps
[`load_go_human_data()`](https://gregorlueg.github.io/bixverse/reference/load_go_human_data.md)
and
[`process_go_data()`](https://gregorlueg.github.io/bixverse/reference/process_go_data.md)
into one.

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

## Examples

``` r
# \donttest{
# human GO data ready for GeneOntologyElim()
go_dt <- get_go_data_human(.verbose = FALSE)
dim(go_dt)
#> [1] 18841     6
# }
```
