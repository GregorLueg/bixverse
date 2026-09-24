# Transfer labels from a Symphony reference to a query via kNN vote

**\[experimental\]** Finds the reference neighbours of every query cell
and assigns the majority label. Ties go to the lowest label index.

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

  Numerical matrix. Reference Harmony-corrected embedding (N_ref x d).

- query_z_corr:

  Numerical matrix. Query Symphony-corrected embedding (N_q x d).

- reference_labels:

  Integer vector. 0-based integer-encoded reference labels.

- n_labels:

  Integer. Number of distinct labels.

- knn_params:

  List. Output of
  [`params_sc_knn()`](https://gregorlueg.github.io/bixverse/reference/params_sc_knn.md).

- seed:

  Integer. Seed for the kNN search.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with

- predicted - Integer vector. Predicted label per query cell (0-based).

- confidence - Numerical vector. Vote share of the winning label.
