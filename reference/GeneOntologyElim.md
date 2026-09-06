# Gene Ontology data

This class is used to store the gene ontology information for usage in
GSE elimination methods.

## Usage

``` r
GeneOntologyElim(go_data_dt, min_genes)
```

## Arguments

- go_data_dt:

  A data.table that contains the gene ontology information. This can be
  extract with for example
  [`get_go_data_human()`](https://gregorlueg.github.io/bixverse/reference/get_go_data_human.md).

- min_genes:

  data.frame. Meta-data information in form of a data.frame.

## Value

Returns the class for subsequent usage.

## Properties

- go_info:

  data.table. Contains the gene ontology identifiers and names.

- go_to_genes:

  List. Contains the genes within each gene ontology term.

- ancestry:

  List. Contains the ancestors for each gene ontology term.

- levels:

  List. Which gene ontology terms sit at which level.

- min_genes:

  Integer, the minimum genes in the gene ontology term to conduct the
  test.

## Examples

``` r
# \donttest{
# human GO restricted to terms with at least 25 genes
go_obj <- GeneOntologyElim(
  get_go_data_human(.verbose = FALSE),
  min_genes = 25L
)
go_obj
#> Gene ontology enrichment class:
#>  Contains 1974 gene ontology terms.
#>  Total of 14 levels represented in the ontology.
#>  Minimum genes per term set to 25.
# }
```
