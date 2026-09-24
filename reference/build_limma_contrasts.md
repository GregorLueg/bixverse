# Build pairwise limma contrasts

Stands in for
[`limma::makeContrasts()`](https://rdrr.io/pkg/limma/man/makeContrasts.html)
for the two cases bixverse needs: every pairwise difference between the
levels of the main contrast, or the differences given as `"a-b"`
strings.

## Usage

``` r
build_limma_contrasts(coef_names, contrast_grps, contrast_list = NULL)
```

## Arguments

- coef_names:

  String vector. The column names of the design matrix.

- contrast_grps:

  String vector. The levels of the main contrast. Only used if
  `contrast_list` is `NULL`.

- contrast_list:

  Optional string vector of the form `"a-b"`. If `NULL`, all pairwise
  contrasts between `contrast_grps` are built, in design column order.

## Value

A named list of numeric contrast vectors, one entry per design column.
Names are the contrasts with `-` replaced by `_vs_`.
