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

## Examples

``` r
# the human gene ontology tables shipped with the package
go_data <- load_go_human_data()
names(go_data)
#> [1] "go_info"       "gene_ontology" "go_to_genes"  
head(go_data$go_info)
#>         go_id                                                  go_name
#>        <char>                                                   <char>
#> 1: GO:0000001                                mitochondrion inheritance
#> 2: GO:0000002                         mitochondrial genome maintenance
#> 3: GO:0000003                                    obsolete reproduction
#> 4: GO:0000005                    obsolete ribosomal chaperone activity
#> 5: GO:0000006    high-affinity zinc transmembrane transporter activity
#> 6: GO:0000007 low-affinity zinc ion transmembrane transporter activity
#>    namespace
#>       <char>
#> 1:        BP
#> 2:        BP
#> 3:        BP
#> 4:        MF
#> 5:        MF
#> 6:        MF
```
