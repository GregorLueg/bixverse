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
