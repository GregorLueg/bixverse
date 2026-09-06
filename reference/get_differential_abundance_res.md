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
#>    Nhood        logFC   logCPM            F    PValue       FDR SpatialFDR
#>    <int>        <num>    <num>        <num>     <num>     <num>      <num>
#> 1:     1  0.696761445 14.32139 8.523319e-01 0.3565560 0.8319639  0.8269572
#> 2:     2  0.340283155 14.32108 2.207650e-01 0.6387624 0.8455142  0.8402161
#> 3:     3 -1.087351007 14.31886 2.594503e+00 0.1267721 0.8163203  0.8108069
#> 4:     4 -0.005885377 14.31889 8.100235e-05 0.9931983 0.9945522  0.9945522
#> 5:     5  0.696752491 14.30611 3.770292e-01 0.5396130 0.8455142  0.8402161
#> 6:     6 -0.351990673 14.32160 2.136750e-01 0.6442013 0.8455142  0.8402161

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
