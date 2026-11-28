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
#>    Nhood       logFC   logCPM            F    PValue       FDR SpatialFDR
#>    <int>       <num>    <num>        <num>     <num>     <num>      <num>
#> 1:     1  0.42390689 13.85440 0.3958957378 0.5296492 0.9455350  0.9459114
#> 2:     2  0.11686896 14.17130 0.0458927514 0.8305024 0.9455350  0.9459114
#> 3:     3 -0.42239717 13.76630 0.3524301245 0.5531450 0.9455350  0.9459114
#> 4:     4 -0.37967351 13.85711 0.3177198405 0.5733617 0.9455350  0.9459114
#> 5:     5  0.72784199 13.61784 1.0455122319 0.3078362 0.9455350  0.9459114
#> 6:     6  0.02071073 13.71601 0.0007936764 0.9775417 0.9775417  0.9775417

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
