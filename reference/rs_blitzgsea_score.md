# Score gene sets against a calibrated blitzGSEA null

**\[experimental\]**

Each gene set costs one enrichment score plus one gamma tail evaluation.
The gene sets are expected to have been filtered to the desired size
bounds and intersected with the signature already.

## Usage

``` r
rs_blitzgsea_score(stats, pathways, null_model, blitz_params)
```

## Arguments

- stats:

  Numeric vector. The gene level statistic. Needs to be sorted in
  descending nature and be the same signature the null was calibrated
  on.

- pathways:

  List. One integer vector of index positions per gene set, indexed to
  R's 1-indexing. Order and duplicates do not matter.

- null_model:

  List. The calibrated null from
  [`rs_blitzgsea_calibrate()`](https://gregorlueg.github.io/bixverse/reference/rs_blitzgsea_calibrate.md).

- blitz_params:

  List. The blitzGSEA parameters, see
  [`params_blitzgsea()`](https://gregorlueg.github.io/bixverse/reference/params_blitzgsea.md).
  Only `centre` is read here and it has to match what the calibration
  used.

## Value

List with the following elements

- es Numeric vector. Enrichment scores for the gene sets.

- nes Numeric vector. Normalised enrichment scores.

- pvals Numeric vector. Two-sided p-values from the gamma approximation.

- sidak Numeric vector. Sidak-adjusted p-values.

- fdr Numeric vector. Benjamini-Hochberg adjusted p-values.

- size Integer vector. Gene set size after intersection.

- leading_edge List of integer vectors with the leading edge index
  positions, indexed to R's 1-indexing.
