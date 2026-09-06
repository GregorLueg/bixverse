# Get the gene names based on the gene idx

Get the gene names based on the gene idx

## Usage

``` r
get_gene_names_from_idx(x, gene_idx, rust_based = TRUE)

# S3 method for class 'ScMap'
get_gene_names_from_idx(x, gene_idx, rust_based = TRUE)

## S7 method for class <bixverse::SingleCells>
get_gene_names_from_idx(x, gene_idx, rust_based = TRUE)

## S7 method for class <bixverse::SingleCellsSubset>
get_gene_names_from_idx(x, gene_idx, rust_based = TRUE)
```

## Arguments

- x:

  An object to get the gene names from.

- gene_idx:

  Integer. The original gene indices for which to return the gene names.

- rust_based:

  Boolean. Is it Rust-based, i.e., 0-index or R-based, i.e., 1-indexed.

## Examples

``` r
# Rust indices translated back into gene identifiers
sc <- demo_single_cells(prepped = FALSE)
get_gene_names_from_idx(sc, gene_idx = 0:2, rust_based = TRUE)
#>         0         1         2 
#> "gene_01" "gene_02" "gene_03" 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
