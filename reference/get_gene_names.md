# Get the gene names

Get the main gene names (for example symbols or Ensembl identifiers).

## Usage

``` r
get_gene_names(x)

# S3 method for class 'ScMap'
get_gene_names(x)

## S7 method for class <bixverse::SingleCells>
get_gene_names(x)

## S7 method for class <bixverse::SingleCellsSubset>
get_gene_names(x)
```

## Arguments

- x:

  An object to get the gene names from.

## Value

The primary gene identifiers stored in the class.

## Examples

``` r
# the primary gene identifiers held by the object
sc <- demo_single_cells(prepped = FALSE)
head(get_gene_names(sc), 3)
#> [1] "gene_01" "gene_02" "gene_03"

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
