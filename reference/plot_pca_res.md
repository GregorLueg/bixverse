# Plot the PCA data

Plot the PCA data

## Usage

``` r
plot_pca_res(
  object,
  cols_to_plot = c("contrast_info", "sample_source"),
  pcs_to_plot = c("PC_1", "PC_2"),
  ...
)
```

## Arguments

- object:

  The underlying class, see [`bulk_dge()`](bulk_dge.md).

- cols_to_plot:

  String vector. The columns within the meta-data to plot. Defaults to
  `c('contrast_info', 'sample_source')`

- pcs_to_plot:

  String vector of length 2. Will default to `c("PC_1", "PC_2")`.

- ...:

  additional parameters

## Value

A plot if the PCA information was found. `NULL` if no PCA was found.
