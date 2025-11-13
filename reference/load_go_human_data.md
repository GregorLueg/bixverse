# Get the Gene Ontology data human

This function loads in gene ontology data stored in the package. This is
for humans only.

## Usage

``` r
load_go_human_data()
```

## Value

A list containing:

- go_info - data.table. The gene ontology identifier, name and namespace
  can be found in this one

- gene_ontology - data.table. The relationships between different gene
  ontology terms.

- go_to_genes - data.table. The gene ontology term to gene (ensembl id)
  relationships.
