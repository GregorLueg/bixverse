# Set gene mapping

Set a gene mapping for a given object. This is used for the single
cell-related classes with streaming from disk.

## Usage

``` r
set_gene_mapping(x, gene_map)

# S3 method for class 'ScMap'
set_gene_mapping(x, gene_map)

## S7 method for class <bixverse::SingleCells>
set_gene_mapping(x, gene_map)
```

## Arguments

- x:

  An object to set gene mapping for

- gene_map:

  Named integer indicating indices and names of the genes

## Examples

``` r
# the mapping is normally written during ingestion
sc <- demo_single_cells(prepped = FALSE)
genes <- get_gene_names(sc)
sc <- set_gene_mapping(sc, stats::setNames(seq_along(genes), genes))
head(get_gene_names(sc), 3)
#> [1] "gene_01" "gene_02" "gene_03"

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
