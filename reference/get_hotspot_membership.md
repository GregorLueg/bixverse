# Get the hotspot gene membership table

Get the hotspot gene membership table

## Usage

``` r
get_hotspot_membership(x)

# S3 method for class 'Hotspot'
get_hotspot_membership(x)
```

## Arguments

- x:

  The object from which to retrieve the hotspot gene membership

## Examples

``` r
# the gene to module table hotspot clustered out
sc <- demo_single_cells()
hs <- hotspot_gene_cor_sc(sc, .verbose = FALSE)
hs <- generate_hotspot_membership(hs)
head(get_hotspot_membership(hs))
#>    gene_id cluster_member
#>     <char>          <num>
#> 1: gene_01              1
#> 2: gene_02              1
#> 3: gene_03              1
#> 4: gene_04              1
#> 5: gene_05              1
#> 6: gene_06              1

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
