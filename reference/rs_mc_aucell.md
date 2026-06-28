# Calculate AUCell in Rust (for meta cells)

**\[experimental\]** The function will take in a list of gene set
indices (0-indexed!) and calculate an AUCell type statistic. Two options
here: calculate this with proper AUROC calculations (useful for marker
gene expression) or based on the Mann-Whitney statistic (useful for
pathway activity measurs). This version works on MetaCell counts which
are stored in memory directly.

## Usage

``` r
rs_mc_aucell(sparse_data, gs_list, auc_type, verbose)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `format`.

- gs_list:

  List. List with the gene set indices (0-indexed!) of the genes of
  interest.

- auc_type:

  String. One of `"wilcox"` or `"auroc"`, pending on which statistic you
  wish to calculate.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A matrix of cells x gene sets with the values representing the AUC.
