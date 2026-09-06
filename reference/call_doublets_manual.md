# Manually readjust Scrublet doublet call thresholds

Updates doublet calls and associated summary statistics in a
`ScrubletRes` object using a user-supplied threshold. Intended for use
after visual inspection of the score histograms via
[`plot.ScrubletRes()`](https://gregorlueg.github.io/bixverse/reference/plot.ScrubletRes.md).

## Usage

``` r
call_doublets_manual(
  scrublet_res,
  threshold,
  for_sample = NULL,
  .verbose = TRUE
)
```

## Arguments

- scrublet_res:

  A `ScrubletRes` object.

- threshold:

  Numeric in `[0, 1]`. The new threshold to apply.

- for_sample:

  Optional character. For grouped results, the name of the group to
  update. Defaults to the first group if `NULL`. Ignored for ungrouped
  results.

- .verbose:

  Logical. If `TRUE`, prints updated rate summaries to the console.

## Value

The `ScrubletRes` object with updated `predicted_doublets`, `z_scores`,
`threshold`, `detected_doublet_rate`, `detectable_doublet_fraction`, and
`overall_doublet_rate`.

## Examples

``` r
# move the automatic threshold after eyeballing the histograms
sc <- demo_single_cells(prepped = FALSE)
res <- scrublet_sc(
  sc,
  scrublet_params = params_scrublet(
    pca = list(no_pcs = 10L),
    hvg = list(min_gene_var_pctl = 0.0),
    n_bins = 20L
  ),
  .verbose = FALSE
)
call_doublets_manual(res, threshold = 0.3, .verbose = FALSE)
#> ScrubletRes: 500 cells, 5 doublets (1.0%)
#>   Threshold:              0.3000
#>   Detected doublet rate:  1.0%
#>   Detectable fraction:    85.3%
#>   Overall doublet rate:   1.2%
#>   Simulated doublets:     750

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
