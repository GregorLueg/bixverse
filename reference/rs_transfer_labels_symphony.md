# Transfer labels from a Symphony reference to a query via kNN vote

**\[experimental\]**

## Usage

``` r
rs_transfer_labels_symphony(
  reference_z_corr,
  query_z_corr,
  reference_labels,
  n_labels,
  knn_params,
  seed,
  verbose
)
```

## Arguments

- reference_z_corr:

  Reference Harmony-corrected embedding (N_ref x d).

- query_z_corr:

  Query Symphony-corrected embedding (N_q x d).

- reference_labels:

  0-based integer-encoded reference labels.

- n_labels:

  Number of distinct labels.

- knn_params:

  List. Output of
  [`params_sc_knn()`](https://gregorlueg.github.io/bixverse/reference/params_sc_knn.md).

- seed:

  Integer.

- verbose:

  Integer. 0/1/2.

## Value

A list with `predicted` (0-based integer per query cell) and
`confidence` (vote share of the winning label).
