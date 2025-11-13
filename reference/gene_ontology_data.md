# Gene Ontology data

This class is used to store the gene ontology information for usage in
GSE elimination methods.

## Usage

``` r
gene_ontology_data(go_data_dt, min_genes)
```

## Arguments

- go_data_dt:

  A data.table that contains the gene ontology information. This can be
  extract with for example
  [`get_go_data_human()`](get_go_data_human.md).

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
