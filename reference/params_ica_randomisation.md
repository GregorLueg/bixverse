# Wrapper function for ICA randomisation

Wrapper function for ICA randomisation

## Usage

``` r
params_ica_randomisation(
  cross_validate = FALSE,
  random_init = 50L,
  folds = 10L
)
```

## Arguments

- cross_validate:

  Boolean. Do you want to apply a cross-validation type approach and
  split the data into `folds` folds to assess within data stability of
  the component.

- random_init:

  Integer. Number of random initialisations to use.

- folds:

  Integer. Number of folds to use if `cross_validate` is set to `TRUE`.
  To note, you will be running `random_init * folds` ICA runs.

## Value

A list with the parameters for usage in the subsequent functions.
