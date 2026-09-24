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
  the component. Defaults to `FALSE`.

- random_init:

  Integer. Number of random initialisations to use. Defaults to `50L`.

- folds:

  Integer. Number of folds to use if `cross_validate` is set to `TRUE`.
  To note, you will be running `random_init * folds` ICA runs. Defaults
  to `10L`.

## Value

A named list with the following elements:

- cross_validate - Boolean. Do you want to apply a cross-validation type
  approach and split the data into `folds` folds to assess within data
  stability of the component. Defaults to `FALSE`.

- random_init - Integer. Number of random initialisations to use.
  Defaults to `50L`.

- folds - Integer. Number of folds to use if `cross_validate` is set to
  `TRUE`. To note, you will be running `random_init * folds` ICA runs.
  Defaults to `10L`.
