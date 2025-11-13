# Process Gene Ontology data into the right format

Helper function that takes in different files containing gene ontology
data and puts them together for various gene set enrichment methods
using the ontological information

## Usage

``` r
process_go_data(go_info, go_genes, go_relationships)
```

## Arguments

- go_info:

  data.table. Contains `go_id`, `go_name` and `namespace.`

- go_genes:

  data.table. Contains `go_id` and corresponding `ensembl_id`.

- go_relationships:

  data.table. Contains `parent`, `child` and `relationship`

## Value

data.table ready for usage in
[`gene_ontology_data()`](gene_ontology_data.md).
