# Calculate Eigengenes for CoReMo modules

This function will calculate the eigengene values for the modules on a
per sample basis and add correlations of the gene expression of a given
gene within the module to its eigengene.

## Usage

``` r
cor_module_coremo_eigengene(object, min_stability = NULL, .verbose = TRUE)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md).

- min_stability:

  Optional float. The minimum stability for the gene you wish to filter
  for based on the leave-one-out resampling. If `NULL`, no filtering
  will be applied.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added correlations to the modules and the values for a
given eigengene per sample as a data.table.

## References

Srivastava, et al., Nat. Commun., 2018; Francois, Romagnolo, et al.,
Nat. Commun., 2024.
