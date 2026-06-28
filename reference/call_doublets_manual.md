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
