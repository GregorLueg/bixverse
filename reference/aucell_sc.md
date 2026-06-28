# Calculate AUC scores (akin to AUCell)

Calculates an AUC-type score akin to AUCell across the gene sets, see
Aibar et al. You have the options to calculate the AUC. Two options
here: calculate this with proper AUROC calculations (useful for marker
gene expression, use the `"auroc"` version) or based on the Mann-Whitney
statistic (useful for pathway activity measurs, use the `"wilcox"`).
Data can be streamed in chunks of 50k cells per or loaded in in one go.

## Usage

``` r
aucell_sc(
  object,
  gs_list,
  auc_type = c("wilcox", "auroc"),
  streaming = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` (or potentially other) class.

- gs_list:

  Named list. The elements have the gene identifiers of the respective
  gene sets.

- auc_type:

  String. Which type of AUC to calculate. Choice of
  `c("wilcox", "auroc")`.

- streaming:

  Optional Boolean. Shall the data be streamed in. Useful for larger
  data sets where you wish to avoid loading in the whole data. If
  `NULL`, will automatically detect. Ignored when applied to
  `MetaCells`.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

AUCell results in form of a matrix that is cells x gene sets or as
`ScMatrixRes` pending the input.

## References

Aibar, et al., Nat Methods, 2017
