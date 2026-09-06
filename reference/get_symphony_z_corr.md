# Getter for the corrected embedding of a Symphony reference

Getter for the corrected embedding of a Symphony reference

## Usage

``` r
get_symphony_z_corr(object)
```

## Arguments

- object:

  `SymphonyReference` class.

## Value

The post-Harmony corrected embedding (N x d).

## Examples

``` r
# the Harmony corrected reference embedding
ref <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 500L,
    n_genes = 50L,
    n_batches = 2L
  )
)
symphony_ref <- build_symphony_ref(
  ref,
  batch_column = "batch_index",
  hvg = get_hvg(ref) + 1L,
  harmony_params = params_sc_harmony(k = 10L),
  no_pcs = 10L,
  label_columns = "cell_grp",
  .verbose = FALSE
)
dim(get_symphony_z_corr(symphony_ref))
#> [1] 500  10

unlink(ref@dir_data, recursive = TRUE, force = TRUE)
```
