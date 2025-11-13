# Run stabilised ICA with a given number of components

This function runs stabilised ICA with the defined number of components.

## Usage

``` r
ica_stabilised_results(
  object,
  no_comp = NULL,
  ica_type = c("logcosh", "exp"),
  iter_params = params_ica_randomisation(),
  ica_params = params_ica_general(),
  random_seed = 42L,
  consistent_sign = TRUE,
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md). You need to apply
  [`ica_processing()`](ica_processing.md) before running this function.

- no_comp:

  Optional integer. Number of components you wish to use for the ICA
  run. If you have run [`ica_evaluate_comp()`](ica_evaluate_comp.md) the
  optimal number is identified via the elbow method and will be used if
  set to `NULL`. You can overwrite this however.

- ica_type:

  String, element of `c("logcosh", "exp")`.

- iter_params:

  List. This list controls the randomisation parameters for the ICA
  runs, see [`params_ica_randomisation()`](params_ica_randomisation.md)
  for estimating stability. Has the following elements:

  - cross_validate - Boolean. Shall the data be split into different
    chunks on which ICA is run. This will slow down the function
    substantially, as every chunk needs to whitened again.

  - random_init - Integer. How many random initialisations shall be used
    for the ICA runs.

  - folds - If `cross_validate` is set to `TRUE` how many chunks shall
    be used. To note, you will run per ncomp random_init \* fold ICA
    runs which can quickly increase.

- ica_params:

  List. The ICA parameters, see
  [`params_ica_general()`](params_ica_general.md) wrapper function. This
  function generates a list containing:

  - maxit - Integer. Maximum number of iterations for ICA.

  - alpha - Float. The alpha parameter for the logcosh version of ICA.
    Should be between 1 to 2.

  - max_tol - Maximum tolerance of the algorithm.

  - verbose - Controls verbosity of the function.

- random_seed:

  Integer. For reproducibility.

- consistent_sign:

  Boolean. If set to `TRUE`, for each source the absolute maximum value
  will be positive, i.e., the sign will be inverted so that the absolute
  bigger tail is set to positive floats.

- .verbose:

  Boolean. Controls verbosity.

## Value

`bulk_coexp` with the the source matrix S, mixing matrix A and other
parameters added to the slots.
