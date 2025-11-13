# Bixverse implementation of ssGSEA

Implementation of the bixverse version of the single sample gene set
enrichment analysis (ssGSEA)

## Usage

``` r
calc_ssgsea(exp, pathways, ssgsea_params = params_ssgsea(), .verbose = FALSE)
```

## Arguments

- exp:

  Numerical matrix. Rows represents the features, columns the
  features/genes.

- pathways:

  List. A named list with each element containing the genes for this
  pathway.

- ssgsea_params:

  List. The GSVA parameters, see [`params_ssgsea()`](params_ssgsea.md)
  wrapper function. This function generates a list containing:

  - alpha - Float. The exponent defining the weight of the tail in the
    random walk performed by ssGSEA.

  - min_size - Integer. Minimum size for the gene sets.

  - max_size - Integer. Maximum size for the gene sets.

  - normalise - Boolean. Shall the scores be normalised.

- .verbose:

  Boolean. Controls verbosity.

## Value

A matrix of shape pathways (that passed the thresholds) x samples.
