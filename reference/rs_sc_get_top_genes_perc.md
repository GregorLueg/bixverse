# Calculates the cumulative proportion of the top X genes

This calculates the cumulative proportion of the top X genes, for
example Top10, 50, 100. High values here indicate low complexity
samples, i.e., bad quality.

## Usage

``` r
rs_sc_get_top_genes_perc(
  f_path_cell,
  top_n_vals,
  cell_indices,
  streaming,
  verbose
)
```

## Arguments

- f_path_cell:

  String. Path to the `counts_cells.bin` file.

- top_n_vals:

  Integer. The different Top X to look for.

- cell_indices:

  Integer. The indices of the cells for which to calculate the
  proportions. (0-indexed!)

- streaming:

  Boolean. Shall the data be worked on in chunks.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the cumulative percentages of the Top X genes defined as in
`top_n_vals`.
