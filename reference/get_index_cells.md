# Get the index cells

Get the index cells

## Usage

``` r
get_index_cells(x)

# S3 method for class 'miloR'
get_index_cells(x)
```

## Arguments

- x:

  An object from which get the index cells.

## Value

The indices of the cells in the neighbourhood.

## Examples

``` r
# the cells the neighbourhoods were centred on
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
head(get_index_cells(milo))
#> [1]  3  8  9 16 28 33

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
