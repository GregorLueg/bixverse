# Calculate AUCell in Rust

**\[experimental\]** The function will take in a list of gene set
indices (0-indexed!) and calculate an AUCell type statistic. Three
options here: the recovery-curve AUC of Aibar, et al. (the actual AUCell
statistic), an AUC derived from the Mann-Whitney statistic, or average
precision. Data can be streamed in chunks of 50k cells per or loaded in
in one go.

## Usage

``` r
rs_aucell(f_path, gs_list, cells_to_keep, aucell_params, streaming, verbose)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- gs_list:

  List. List with the gene set indices (0-indexed!) of the genes of
  interest.

- cells_to_keep:

  Integer. Vector of indices of the cells to keep.

- aucell_params:

  List. The AUCell parameters, see
  [`params_sc_aucell()`](https://gregorlueg.github.io/bixverse/reference/params_sc_aucell.md).

- streaming:

  Boolean. Shall the data be streamed.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A matrix of cells x gene sets with the values representing the AUC.
