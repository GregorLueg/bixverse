# Split CoReMo modules by correlation sign

This function will split the identified modules by their correlation
sign. In certain cases, positive and negatively correlated genes can be
part of the same module. This function will annotate them with `"_pos"`
and `"_neg"` respectively.

## Usage

``` r
cor_module_coremo_cor_sign(object, min_stability = NULL, .verbose = TRUE)
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

The class with updated correlation module names.

## References

Srivastava, et al., Nat. Commun., 2018; Francois, Romagnolo, et al.,
Nat. Commun., 2024.
