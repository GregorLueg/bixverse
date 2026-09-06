# Getter for the stored labels of a Symphony reference

Getter for the stored labels of a Symphony reference

## Usage

``` r
get_symphony_labels(object)
```

## Arguments

- object:

  `SymphonyReference` class.

## Value

A `data.table` of reference cell labels in `z_corr` row order, or `NULL`
if no labels are stored.

## Examples

``` r
# the cell labels snapshotted at build time
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
head(get_symphony_labels(symphony_ref))
#>       cell_grp
#>         <char>
#> 1: cell_type_1
#> 2: cell_type_2
#> 3: cell_type_3
#> 4: cell_type_1
#> 5: cell_type_2
#> 6: cell_type_3

unlink(ref@dir_data, recursive = TRUE, force = TRUE)
```
