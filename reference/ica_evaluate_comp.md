# Iterate over different ncomp parameters for ICA

This function allows to iterate over a vector of ncomp to identify which
ncomp parameter to choose for your data set. The idea is to generate
stability profiles over the different ncomps and identify a 'sweet spot'
of good stability, low mutual information and good convergence of the
identified independent components

## Usage

``` r
ica_evaluate_comp(
  object,
  ica_type = c("logcosh", "exp"),
  iter_params = params_ica_randomisation(),
  ncomp_params = params_ica_ncomp(),
  ica_params = params_ica_general(),
  random_seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md). You need to apply
  [`ica_processing()`](ica_processing.md) before running this function.

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

- ncomp_params:

  List. Parameters for the ncomp to iterate through, see
  [`params_ica_ncomp()`](params_ica_ncomp.md). In the standard setting,
  `c(2, 3, 4, 5)` will be tested and then in steps until max_no_comp
  will be tested, i.e.,
  `c(2, 3, 4, 5, 10, 15, ..., max_no_comp - 5, max_no_comp)`.

  - max_no_comp - Maximum number of ncomp to test.

  - steps - Integer. In which steps to move from 5 onwards.

  - custom_seq - An integer vector. If you wish to provide a custom
    version of no_comp to iterate through.

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

- .verbose:

  Boolean. Controls verbosity.

## Value

`bulk_coexp` with the added information of stability of the components
and other data to plot to choose the right `ncomp`.
