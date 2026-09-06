# Getter the var table

Getter the var table

## Usage

``` r
get_sc_var(object, indices = NULL, cols = NULL, modality = c("rna", "adt"))
```

## Arguments

- object:

  `SingleCells`, `MetaCells`, `SingleCellsMultiModal` class.

- indices:

  Optional integer vector. The integer positions of the genes to return.

- cols:

  Optional string vector. The columns from the var table to return.

- modality:

  String. The modality to return. One of `c("rna", "adt")`.

## Value

The vars table

## Examples

``` r
# the per gene statistics the HVG step wrote
sc <- demo_single_cells()
head(get_sc_var(sc, cols = c("gene_id", "mean", "var_std")), 3)
#>    gene_id   mean   var_std
#>     <char>  <num>     <num>
#> 1: gene_01 11.922 0.9063808
#> 2: gene_02 10.294 0.8359538
#> 3: gene_03 12.002 0.7752684

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
