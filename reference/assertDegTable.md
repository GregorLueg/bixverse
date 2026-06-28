# Assert a DE table has the expected columns and types

Assert a DE table has the expected columns and types

## Usage

``` r
assertDegTable(x, name, needs_cluster = TRUE)
```

## Arguments

- x:

  The candidate `data.table`.

- name:

  Variable name for error messages.

- needs_cluster:

  If TRUE, also require a `cluster_id` column.
