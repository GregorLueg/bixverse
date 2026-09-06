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
#>    Nhood        logFC   logCPM            F    PValue       FDR SpatialFDR
#>    <int>        <num>    <num>        <num>     <num>     <num>      <num>
#> 1:     1  0.696761445 14.32139 8.523319e-01 0.3565560 0.8319639  0.8269572
#> 2:     2  0.340283155 14.32108 2.207650e-01 0.6387624 0.8455142  0.8402161
#> 3:     3 -1.087351007 14.31886 2.594503e+00 0.1267721 0.8163203  0.8108069
#> 4:     4 -0.005885377 14.31889 8.100235e-05 0.9931983 0.9945522  0.9945522
#> 5:     5  0.696752491 14.30611 3.770292e-01 0.5396130 0.8455142  0.8402161
#> 6:     6 -0.351990673 14.32160 2.136750e-01 0.6442013 0.8455142  0.8402161
#>    majority_celltype majority_prop
#>               <char>         <num>
#> 1:       cell_type_1        1.0000
#> 2:       cell_type_3        1.0000
#> 3:       cell_type_1        0.6875
#> 4:       cell_type_2        0.9375
#> 5:       cell_type_2        1.0000
#> 6:       cell_type_1        1.0000

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
