# Identify hotspot gene clusters

Identify hotspot gene clusters

## Usage

``` r
generate_hotspot_membership(x, fdr_threshold = 0.05, min_size = 10L)

# S3 method for class 'Hotspot'
generate_hotspot_membership(x, fdr_threshold = 0.05, min_size = 10L)
```

## Arguments

- x:

  An object to generate the hotspot gene clusters for.

- fdr_threshold:

  Numeric. The maximum FDR for a given gene-gene local correlation to be
  included.

- min_size:

  Integer. Minimum cluster size.

## Examples

``` r
# cluster the local gene-gene correlations into modules
sc <- demo_single_cells()
hs <- hotspot_gene_cor_sc(sc, .verbose = FALSE)
generate_hotspot_membership(hs)
#> Hotspot gene-gene local correlation results
#>   Genes: 50
#>   Cells: 500
#>   Modules: 3 (31 genes assigned, 19 unassigned)

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
