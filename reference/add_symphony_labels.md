# Add labels to a Symphony reference post-hoc

Reads one or more obs columns from a `SingleCells` object and stores
them in the reference's `labels` slot for use by
[`transfer_labels_symphony()`](https://gregorlueg.github.io/bixverse/reference/transfer_labels_symphony.md).
The provided `sc_object` must have `cells_to_keep` matching the cells
used at reference-build time; this is enforced by a length check against
`nrow(z_corr)`.

## Usage

``` r
add_symphony_labels(reference, sc_object, columns, overwrite = FALSE)
```

## Arguments

- reference:

  `SymphonyReference`.

- sc_object:

  `SingleCells` to read labels from. Typically the same object passed to
  [`build_symphony_ref()`](https://gregorlueg.github.io/bixverse/reference/build_symphony_ref.md).

- columns:

  Character vector of obs column names.

- overwrite:

  Boolean. If `TRUE`, existing label columns of the same name are
  replaced; otherwise an error is raised on collision.

## Value

The `reference` with updated `labels`.

## Examples

``` r
# attach obs labels to a reference built without them
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
  .verbose = FALSE
)
symphony_ref <- add_symphony_labels(
  symphony_ref,
  sc_object = ref,
  columns = "cell_grp"
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
