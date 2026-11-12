# Get the HVG

Returns the positions of the HVGs in the var table. These are 0-based on
every class, as Rust wants them: add 1 before indexing anything in R,
e.g. `get_gene_names(sc)[get_hvg(sc) + 1L]`.

## Usage

``` r
get_hvg(x)

## S7 method for class <bixverse::MetaCells>
get_hvg(x)

# S3 method for class 'ScMap'
get_hvg(x)

## S7 method for class <bixverse::SingleCells>
get_hvg(x)

## S7 method for class <bixverse::SingleCellsSubset>
get_hvg(x)
```

## Arguments

- x:

  An object to get HVG from.

## Value

Integer vector of 0-based gene indices, empty with a warning when no
HVGs are stored.

## Examples

``` r
# stored 0-based, so map them back through the gene names
sc <- demo_single_cells()
get_gene_names_from_idx(sc, head(get_hvg(sc), 3))
#>        43        18         4 
#> "gene_44" "gene_19" "gene_05" 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
