# Calculate AUC scores (akin to AUCell)

Calculates an AUC-type score akin to AUCell across the gene sets. You
have the options to calculate the AUC. Two options here: calculate this
with proper AUROC calculations (useful for marker gene expression, use
the `"auroc"` version) or based on the Mann-Whitney statistic (useful
for pathway activity measurs, use the `"wilcox"`). Data can be streamed
in chunks of 50k cells per or loaded in in one go.

## Usage

``` r
aucell_sc(
  object,
  gs_list,
  auc_type = c("wilcox", "auroc"),
  streaming = FALSE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- gs_list:

  Named list. The elements have the gene identifiers of the respective
  gene sets.

- auc_type:

  String. Which type of AUC to calculate. Choice of
  `c("wilcox", "auroc")`.

- streaming:

  Boolean. Shall the cell data be streamed in. Useful for larger data
  sets.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

data.table with the DGE results from the test.
