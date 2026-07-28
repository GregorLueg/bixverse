# Wrapper function for parameters for AUCell

The three statistics consume the same within-cell ranking, but weight it
very differently. `"wilcox"` is a pure function of the gene set's rank
sum, so a gene at rank 2 and one at rank 200 count for almost the same
thing. `"recovery"` and `"ap"` are top-heavy.

## Usage

``` r
params_sc_aucell(
  auc_type = c("wilcox", "recovery", "ap"),
  max_rank = NULL,
  standardise = FALSE
)
```

## Arguments

- auc_type:

  String. Which statistic to calculate. One of
  `c("wilcox", "recovery", "ap")`. `"wilcox"` is the AUC derived from
  the Mann-Whitney U statistic over the full ranking, with the null at
  0.5 for any gene set size. `"recovery"` is the recovery-curve AUC
  under a rank cutoff, i.e. the actual AUCell statistic of Aibar, et al.
  `"ap"` is average precision, the most top-heavy of the three, but its
  null tracks the gene set prevalence so raw values are not comparable
  across gene sets of different size unless `standardise` is on.
  Defaults to `"wilcox"`.

- max_rank:

  Optional numeric. Rank cutoff for `"recovery"`, counted from the top
  of the within-cell ranking. If `NULL`, resolves to the top 5% of the
  gene universe, following Aibar, et al. Ignored by the other two
  statistics.

- standardise:

  Boolean. Shall each gene set's scores be z-scored across the cells.
  This is what makes `"ap"` comparable across gene sets of different
  size. Defaults to `FALSE`.

## Value

A list with the AUCell parameters.

## References

Aibar, et al., Nat Methods, 2017
