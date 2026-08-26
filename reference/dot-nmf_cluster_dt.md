# Assemble the per-component consensus diagnostics

Turns the flat vectors the Rust side returns into one row per pooled
component. Rust pools component `j` of restart `r` at position
`k * r + j`, which is the order
[`new_stabilised_nmf_result()`](https://gregorlueg.github.io/bixverse/reference/new_stabilised_nmf_result.md)
builds its `w_all` column names in, so `component_id` joins straight
onto the columns of
[`get_w()`](https://gregorlueg.github.io/bixverse/reference/get_w.md) on
a stabilised fit.

## Usage

``` r
.nmf_cluster_dt(nmf_res, k, n_runs)
```

## Arguments

- nmf_res:

  Named list. Raw consensus output, needs `labels`, `local_density`,
  `kept` and `silhouette`.

- k:

  Integer. Components per restart.

- n_runs:

  Integer. Number of restarts.

## Value

A data.table with one row per pooled component and the columns
`component_id`, `run`, `component`, `pooled_idx`, `cluster`,
`local_density`, `silhouette` and `kept`. Dropped components carry `NA`
in `cluster` and `silhouette`.
