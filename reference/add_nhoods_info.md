# Add neighbourhood info on majority cell type

This function adds cell type composition information to the
`nhoods_info` slot within the `miloR` object. For each neighbourhood, it
calculates the proportion of the majority cell type and identifies which
cell type is most abundant. This is useful for annotating differential
abundance results with the cellular composition of each neighbourhood.

## Usage

``` r
add_nhoods_info(x, cell_info)

# S3 method for class 'miloR'
add_nhoods_info(x, cell_info)
```

## Arguments

- x:

  `miloR` object on which to tag on additional neighbourhood
  information.

- cell_info:

  Character vector. Represents the cell type annotations you wish to add
  to the different neighbourhoods. Must be the same length as the number
  of cells (rows) in the nhoods matrix.

## Value

Modified `miloR` object with updated `nhoods_info` containing
`majority_celltype` and `majority_prop` columns.

## Examples

``` r
# tag each neighbourhood with its majority cell type
sc <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 500L,
    n_genes = 50L,
    n_samples = 6L,
    sample_bias = "even"
  )
)
milo <- get_miloR_abundances_sc(
  sc,
  sample_id_col = "sample_id",
  miloR_params = params_sc_miloR(k_refine = 10L),
  .verbose = FALSE
)
design_df <- data.frame(
  grp = rep(c("a", "b"), each = 3),
  row.names = sprintf("sample_%i", 1:6)
)
milo <- test_nhoods(milo, design = ~grp, design_df = design_df)
milo <- add_nhoods_info(milo, cell_info = get_sc_obs(sc)$cell_grp)
head(get_differential_abundance_res(milo))
#>    Nhood       logFC   logCPM            F    PValue       FDR SpatialFDR
#>    <int>       <num>    <num>        <num>     <num>     <num>      <num>
#> 1:     1  0.42390689 13.85440 0.3958957378 0.5296492 0.9455350  0.9459114
#> 2:     2  0.11686896 14.17130 0.0458927514 0.8305024 0.9455350  0.9459114
#> 3:     3 -0.42239717 13.76630 0.3524301245 0.5531450 0.9455350  0.9459114
#> 4:     4 -0.37967351 13.85711 0.3177198405 0.5733617 0.9455350  0.9459114
#> 5:     5  0.72784199 13.61784 1.0455122319 0.3078362 0.9455350  0.9459114
#> 6:     6  0.02071073 13.71601 0.0007936764 0.9775417 0.9775417  0.9775417
#>    majority_celltype majority_prop
#>               <char>         <num>
#> 1:       cell_type_1     1.0000000
#> 2:       cell_type_3     1.0000000
#> 3:       cell_type_1     0.6315789
#> 4:       cell_type_2     0.8571429
#> 5:       cell_type_2     1.0000000
#> 6:       cell_type_1     1.0000000

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
