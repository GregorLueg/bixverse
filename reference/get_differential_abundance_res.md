# Get the differential abundance results

Get the differential abundance results

## Usage

``` r
get_differential_abundance_res(x)

# S3 method for class 'miloR'
get_differential_abundance_res(x)
```

## Arguments

- x:

  An object from which to get the differential abundance results from.

## Value

The differential abundance results stored in the object if found.

## Examples

``` r
# neighbourhood level differential abundance table
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
head(get_differential_abundance_res(milo))
#>    Nhood       logFC   logCPM           F    PValue       FDR SpatialFDR
#>    <int>       <num>    <num>       <num>     <num>     <num>      <num>
#> 1:     1  0.43394889 13.84157 0.433206314 0.5108705 0.9640042  0.9599185
#> 2:     2  0.12629598 14.15793 0.053925002 0.8165104 0.9640042  0.9644670
#> 3:     3 -0.53221669 13.79775 0.599966995 0.4391353 0.9640042  0.9598125
#> 4:     4 -0.36980693 13.84297 0.313858391 0.5756954 0.9640042  0.9599185
#> 5:     5  0.73814666 13.60314 1.138735255 0.2936231 0.9640042  0.9598125
#> 6:     6  0.03387627 13.70417 0.002060917 0.9638176 0.9754298  0.9772304

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
