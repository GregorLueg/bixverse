# Calculate AUCell in Rust (for meta cells)

The function will take in a list of gene set indices (0-indexed!) and
calculate an AUCell type statistic. Two options here: calculate this
with proper AUROC calculations (useful for marker gene expression) or
based on the Mann-Whitney statistic (useful for pathway activity
measurs). This version works on MetaCell counts which are stored in
memory directly.

## Usage

``` r
rs_mc_aucell(sparse_data, gs_list, auc_type, verbose)
```

## Arguments

- gs_list:

  List. List with the gene set indices (0-indexed!) of the genes of
  interest.

- auc_type:

  String. One of `"wilcox"` or `"auroc"`, pending on which statistic you
  wish to calculate.

- verbose:

  Boolean. Controls verbosity of the function.

- cells_to_keep:

  Integer. Vector of indices of the cells to keep.

- streaming:

  Boolean. Shall the data be streamed.

## Value

A matrix of cells x gene sets with the values representing the AUC.
