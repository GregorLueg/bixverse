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
#>    Nhood       logFC   logCPM           F    PValue       FDR SpatialFDR
#>    <int>       <num>    <num>       <num>     <num>     <num>      <num>
#> 1:     1  0.43394889 13.84157 0.433206314 0.5108705 0.9640042  0.9599185
#> 2:     2  0.12629598 14.15793 0.053925002 0.8165104 0.9640042  0.9644670
#> 3:     3 -0.53221669 13.79775 0.599966995 0.4391353 0.9640042  0.9598125
#> 4:     4 -0.36980693 13.84297 0.313858391 0.5756954 0.9640042  0.9599185
#> 5:     5  0.73814666 13.60314 1.138735255 0.2936231 0.9640042  0.9598125
#> 6:     6  0.03387627 13.70417 0.002060917 0.9638176 0.9754298  0.9772304
#>    majority_celltype majority_prop
#>               <char>         <num>
#> 1:       cell_type_1     1.0000000
#> 2:       cell_type_3     1.0000000
#> 3:       cell_type_1     0.6000000
#> 4:       cell_type_2     0.8571429
#> 5:       cell_type_2     1.0000000
#> 6:       cell_type_1     1.0000000

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
