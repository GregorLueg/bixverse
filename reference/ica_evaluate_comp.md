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

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).
  You need to apply
  [`ica_processing()`](https://gregorlueg.github.io/bixverse/reference/ica_processing.md)
  before running this function.

- ica_type:

  String, element of `c("logcosh", "exp")`.

- iter_params:

  List. This list controls the randomisation parameters for the ICA
  runs, see
  [`params_ica_randomisation()`](https://gregorlueg.github.io/bixverse/reference/params_ica_randomisation.md)
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
  [`params_ica_ncomp()`](https://gregorlueg.github.io/bixverse/reference/params_ica_ncomp.md).
  In the standard setting, `c(2, 3, 4, 5)` will be tested and then in
  steps until max_no_comp will be tested, i.e.,
  `c(2, 3, 4, 5, 10, 15, ..., max_no_comp - 5, max_no_comp)`.

  - max_no_comp - Maximum number of ncomp to test.

  - steps - Integer. In which steps to move from 5 onwards.

  - custom_seq - An integer vector. If you wish to provide a custom
    version of no_comp to iterate through.

- ica_params:

  List. The ICA parameters, see
  [`params_ica_general()`](https://gregorlueg.github.io/bixverse/reference/params_ica_general.md)
  wrapper function. This function generates a list containing:

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

`BulkCoExp` with the added information of stability of the components
and other data to plot to choose the right `ncomp`.

## Examples

``` r
# stability across a small grid of component counts
mat <- t(synthetic_signal_matrix()$mat)
obj <- BulkCoExp(mat, data.table::data.table(sample_id = rownames(mat)))
obj <- preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
obj <- ica_processing(obj, .verbose = FALSE)
obj <- ica_evaluate_comp(
  obj,
  ica_type = "logcosh",
  ncomp_params = params_ica_ncomp(custom_seq = seq(2L, 20L, by = 2L)),
  .verbose = FALSE
)
head(get_ica_stability_res(obj))
#> Key: <no_components>
#>    no_components median_stability UQ_stability LQ_stability converged
#>            <int>            <num>        <num>        <num>     <num>
#> 1:             2        0.3341118    0.3371207    0.3311030      1.00
#> 2:             4        0.5896017    0.6422570    0.5311256      1.00
#> 3:             6        0.6287036    0.6477328    0.5436103      0.34
#> 4:             8        0.4903164    0.5704741    0.4185860      0.00
#> 5:            10        0.6038098    0.7605020    0.4277218      0.58
#> 6:            12        0.4654485    0.5977709    0.3903346      0.30
#>    norm_mutual_information combined_score
#>                      <num>          <num>
#> 1:              0.28908019      0.2375267
#> 2:              0.15172833      0.5001424
#> 3:              0.10588218      0.1911259
#> 4:              0.09540071      0.0000000
#> 5:              0.09111733      0.3182995
#> 6:              0.09419421      0.1264818
```
