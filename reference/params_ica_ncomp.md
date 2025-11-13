# Wrapper function for ICA ncomp iterations

Wrapper function to provide parameters through which ncomps to iterate
through.

## Usage

``` r
params_ica_ncomp(max_no_comp = 75L, steps = 5L, custom_seq = NULL)
```

## Arguments

- max_no_comp:

  Integer. Maximum number of ncomp to test.

- steps:

  Integer. In which steps to move from 5 onwards.

- custom_seq:

  An integer vector. If you wish to provide a custom version of no_comp
  to iterate through. If NULL, you will iterate through
  `c(2, 3, 4, 5, 5 + step, ... max_no_comp - step, max_no_comp)`

## Value

A list with the parameters for usage in subsequent functions.
