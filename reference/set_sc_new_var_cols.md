# Add a new column to the var table

Add a new column to the var table

## Usage

``` r
set_sc_new_var_cols(object, data_list)
```

## Arguments

- object:

  `SingleCells` class.

- data_list:

  Named list with the data to add.

## Value

The class with updated var table in the DuckDB

## Examples

``` r
# flag a couple of genes in the var table
sc <- demo_single_cells(prepped = FALSE)
sc <- set_sc_new_var_cols(
  sc,
  data_list = list(
    is_marker = get_gene_names(sc) %in% c("gene_01", "gene_02")
  )
)
head(get_sc_var(sc), 3)
#>    gene_idx gene_id ensembl_id no_cells_exp is_marker
#>       <int>  <char>     <char>        <int>    <lgcl>
#> 1:        1 gene_01     ens_01          408      TRUE
#> 2:        2 gene_02     ens_02          403      TRUE
#> 3:        3 gene_03     ens_03          409     FALSE

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
