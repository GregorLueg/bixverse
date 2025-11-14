# Calculate VISION pathway scores in Rust

The function will take in a list of gene sets that contains lists of
`"pos"` and `"neg"` gene indices (0-indexed). You don't have to provide
the `"neg"`, but it can be useful to classify the delta of two stats
(EMT, Th1; Th2) etc.

## Usage

``` r
rs_vision(f_path, gs_list, cells_to_keep, streaming, verbose)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- gs_list:

  Nested list. Each sublist contains the (0-indexed!) positive and
  negative gene indices of that specific gene set.

- cells_to_keep:

  Integer. Vector of indices of the cells to keep.

- streaming:

  Boolean. Shall the data be streamed.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A matrix of cells x vision scores per gene set.
