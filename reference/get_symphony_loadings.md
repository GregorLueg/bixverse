# Getter for the PCA loadings of a Symphony reference

Getter for the PCA loadings of a Symphony reference

## Usage

``` r
get_symphony_loadings(object)
```

## Arguments

- object:

  `SymphonyReference` class.

## Value

The PCA gene loadings matrix (n_hvgs x d).

## Examples

``` r
# the PCA gene loadings the query gets projected through
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
dim(get_symphony_loadings(symphony_ref))
#> [1] 30 10

unlink(ref@dir_data, recursive = TRUE, force = TRUE)
```
