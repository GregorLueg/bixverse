# Getter for the HVG gene names of a Symphony reference

Getter for the HVG gene names of a Symphony reference

## Usage

``` r
get_symphony_hvg_names(object)
```

## Arguments

- object:

  `SymphonyReference` class.

## Value

Character vector of HVG gene names in reference loading order.

## Examples

``` r
# the HVG names in reference loading order
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
head(get_symphony_hvg_names(symphony_ref))
#> [1] "gene_41" "gene_32" "gene_37" "gene_50" "gene_46" "gene_47"

unlink(ref@dir_data, recursive = TRUE, force = TRUE)
```
