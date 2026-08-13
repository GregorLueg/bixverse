# Calculate AUC scores (akin to AUCell)

Calculates an AUC-type score akin to AUCell across the gene sets, see
Aibar et al. Three statistics are on offer, all consuming the same
within-cell ranking but weighting it differently. `"wilcox"` (the
default) is the AUC derived from the Mann-Whitney U statistic over the
full ranking; its null sits at 0.5 for any gene set size, which makes it
a good fit for pathway activity. `"recovery"` is the recovery-curve AUC
under a rank cutoff, i.e. the actual AUCell statistic of Aibar, et al.,
and is top-heavy: only genes inside the top `max_rank` of the cell
contribute. `"ap"` is average precision, the most top-heavy of the
three, but its null tracks the gene set prevalence, so raw values are
not comparable across gene sets of different size unless `standardise`
is on. Data can be streamed in chunks of 50k cells per or loaded in in
one go.

## Usage

``` r
aucell_sc(
  object,
  gs_list,
  aucell_params = params_sc_aucell(),
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

- aucell_params:

  List with the AUCell parameters, see
  [`params_sc_aucell()`](https://gregorlueg.github.io/bixverse/reference/params_sc_aucell.md)
  with the following elements:

  - auc_type - String. Which statistic to calculate. One of
    `c("recovery", "wilcox", "ap")`. `"recovery"` is the SCENIC one.

  - max_rank - Optional numeric. Rank cutoff for `"recovery"`. If
    `NULL`, the top 5% of the gene universe is used. Ignored by the
    other statistics.

  - standardise - Boolean. Shall each gene set's scores be z-scored
    across the cells.

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
