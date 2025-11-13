# Create all limma contrasts (based on combination of everything)

Create all limma contrasts (based on combination of everything)

## Usage

``` r
all_limma_contrasts(limma_fit, contrast_grps)
```

## Arguments

- limma_fit:

  The fitted limma model, i.e., output of
  [`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html).

- contrast_grps:

  String vector. The contrast groups of interest. If NULL all co-variate
  comparisons will be returned.

## Value

The Limma contrasts for further usage.
