# Gene Ontology data (deprecated)

**\[deprecated\]**

This constructor has been renamed to
[`GeneOntologyElim()`](https://gregorlueg.github.io/bixverse/reference/GeneOntologyElim.md).

## Usage

``` r
gene_ontology_data(go_data_dt, min_genes)
```

## Arguments

- go_data_dt:

  A data.table that contains the gene ontology information. This can be
  extracted with for example
  [`get_go_data_human()`](https://gregorlueg.github.io/bixverse/reference/get_go_data_human.md).

- min_genes:

  Integer. The minimum number of genes in a gene ontology term to
  conduct the test.

## Value

Returns a
[`GeneOntologyElim()`](https://gregorlueg.github.io/bixverse/reference/GeneOntologyElim.md)
object.
