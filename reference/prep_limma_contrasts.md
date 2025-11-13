# Create all limma contrasts from a provided string

Create all limma contrasts from a provided string

## Usage

``` r
prep_limma_contrasts(limma_fit, contrast_list)
```

## Arguments

- limma_fit:

  The fitted limma model, i.e., output of
  [`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html).

- contrast_list:

  String vector. The strings need to have the form of
  `"contrast1-contrast2"`.

## Value

The Limma contrasts for further usage.
