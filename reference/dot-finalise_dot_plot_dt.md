# Finalise a long dot-plot data.table

Finalise a long dot-plot data.table

## Usage

``` r
.finalise_dot_plot_dt(plot_dt, features, group_levels, scale_exp)
```

## Arguments

- plot_dt:

  data.table. With columns `gene`, `group`, `mean_exp`, `pct_exp`.

- features:

  Character vector. Feature ids, in display order.

- group_levels:

  Character vector. Group levels, in display order.

- scale_exp:

  Boolean. Whether to min-max scale mean expression per gene.

## Value

The data.table with an added `scaled_exp` column and ordered
`gene`/`group` factors.
