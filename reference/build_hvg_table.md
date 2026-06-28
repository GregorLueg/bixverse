# Assemble the HVG data.table and flag the top N

Internal helper shared by the `get_hvg_data_sc` methods. Takes the
per-gene statistics returned by the Rust HVG functions, attaches them to
a copy of the var table, and flags the top `hvg_no` genes according to
the ranking column appropriate for the chosen HVG method.

## Usage

``` r
build_hvg_table(var_table, res, hvg_no, hvg_method)
```

## Arguments

- var_table:

  data.table. The var table (or a subset of it) containing at least
  `gene_idx` and `gene_id`. Copied internally so the caller's table is
  not mutated.

- res:

  Named list. Per-gene statistics returned by `rs_sc_hvg` or
  `rs_mc_hvg`. The names of `res` become new columns on the returned
  data.table. Must contain the ranking column implied by `hvg_method`
  (`var_std` for `"vst"`, `dispersion` for `"dispersion"`,
  `dispersion_scaled` for `"meanvarbin"`).

- hvg_no:

  Integer. Number of top genes to flag as HVGs.

- hvg_method:

  String. One of `c("vst", "dispersion", "meanvarbin")`. Selects which
  column in `res` is used to rank genes.

## Value

A data.table with the original `var_table` columns, all columns from
`res`, plus:

- `is_hvg` - Boolean. `TRUE` for the top `hvg_no` genes.

- `hvg_rank` - Integer. Rank from 1 (most variable) to `hvg_no` for
  HVGs, `NA_integer_` otherwise.
