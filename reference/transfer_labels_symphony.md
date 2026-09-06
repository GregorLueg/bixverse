# Transfer labels from a Symphony reference to a query via kNN majority vote

kNN majority vote on the Harmony-corrected embeddings. The reference's
`z_corr` is the searchable index; the query's `"symphony"` embedding is
what gets queried. Reference labels come from the reference's stored
`labels` slot — populate it via the `label_columns` argument to
[`build_symphony_ref()`](https://gregorlueg.github.io/bixverse/reference/build_symphony_ref.md)
or post-hoc via
[`add_symphony_labels()`](https://gregorlueg.github.io/bixverse/reference/add_symphony_labels.md).

Distances on `z_corr` are typically Euclidean — set
`knn_params$ann_dist` to `"euclidean"` unless you have a reason to use
cosine.

## Usage

``` r
transfer_labels_symphony(
  reference,
  query,
  label_column,
  knn_params = params_sc_knn(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- reference:

  `SymphonyReference` (must have `z_corr` and stored labels).

- query:

  `SingleCells` with a `"symphony"` embedding attached.

- label_column:

  String. Name of a column in the reference's stored labels.

- knn_params:

  List. Output of
  [`params_sc_knn()`](https://gregorlueg.github.io/bixverse/reference/params_sc_knn.md).

- seed:

  Integer.

- .verbose:

  Boolean or integer.

## Value

A data.table with columns `predicted_<label_column>` and
`confidence_<label_column>`, in `get_cells_to_keep(query)` order.

## Examples

``` r
# kNN label transfer from the reference onto a mapped query
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
query <- demo_single_cells(prepped = FALSE, seed = 7L)
query <- map_symphony_query(symphony_ref, query = query, .verbose = FALSE)
labels <- transfer_labels_symphony(
  symphony_ref,
  query = query,
  label_column = "cell_grp",
  .verbose = FALSE
)
head(labels)
#>    predicted_cell_grp confidence_cell_grp
#>                <char>               <num>
#> 1:        cell_type_1           0.8666667
#> 2:        cell_type_2           0.9333333
#> 3:        cell_type_3           0.8000000
#> 4:        cell_type_1           0.8000000
#> 5:        cell_type_2           0.8666667
#> 6:        cell_type_3           0.8666667

unlink(c(ref@dir_data, query@dir_data), recursive = TRUE, force = TRUE)
```
