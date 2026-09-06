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
[`GeneOntologyElim()`](https://gregorlueg.github.io/bixverse/reference/GeneOntologyElim.md).

## Examples

``` r
# \donttest{
# assemble the packaged human GO data by hand
go_data <- load_go_human_data()
relationships <- data.table::setnames(
  data.table::copy(go_data$gene_ontology),
  old = c("from", "to"),
  new = c("parent", "child")
)
go_dt <- process_go_data(
  go_info = go_data$go_info,
  go_genes = go_data$go_to_genes,
  go_relationships = relationships[relationship %in% c("is_a", "part_of")]
)
dim(go_dt)
#> [1] 18841     6
# }
```
