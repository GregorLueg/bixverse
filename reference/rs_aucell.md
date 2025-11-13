# Calculate AUCell in Rust

The function will take in a list of gene set indices (0-indexed!) and
calculate an AUCell type statistic. Two options here: calculate this
with proper AUROC calculations (useful for marker gene expression) or
based on the Mann-Whitney statistic (useful for pathway activity
measurs). Data can be streamed in chunks of 50k cells per or loaded in
in one go.

## Usage

``` r
rs_aucell(f_path, gs_list, cells_to_keep, auc_type, streaming, verbose)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- gs_list:

  List. List with the gene set indices (0-indexed!) of the genes of
  interest.

- cells_to_keep:

  Integer. Vector of indices of the cells to keep.

- auc_type:

  String. One of `"wilcox"` or `"auroc"`, pending on which statistic you
  wish to calculate.

- streaming:

  Boolean. Shall the data be streamed.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A matrix of gene set AUCs x cells.
