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
