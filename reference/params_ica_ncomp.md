# Wrapper function for ICA ncomp iterations

Wrapper function to provide parameters through which ncomps to iterate
through.

## Usage

``` r
params_ica_ncomp(max_no_comp = 75L, steps = 5L, custom_seq = NULL)
```

## Arguments

- max_no_comp:

  Integer. Maximum number of ncomp to test. Defaults to `75L`.

- steps:

  Integer. In which steps to move from 5 onwards. Defaults to `5L`.

- custom_seq:

  Integer vector or `NULL`. If you wish to provide a custom version of
  no_comp to iterate through. If NULL, you will iterate through
  `c(2, 3, 4, 5, 5 + step, ... max_no_comp - step, max_no_comp)`
  Defaults to `NULL`.

## Value

A named list with the following elements:

- max_no_comp - Integer. Maximum number of ncomp to test. Defaults to
  `75L`.

- steps - Integer. In which steps to move from 5 onwards. Defaults to
  `5L`.

- custom_seq - Integer vector or `NULL`. If you wish to provide a custom
  version of no_comp to iterate through. If NULL, you will iterate
  through `c(2, 3, 4, 5, 5 + step, ... max_no_comp - step, max_no_comp)`
  Defaults to `NULL`.
