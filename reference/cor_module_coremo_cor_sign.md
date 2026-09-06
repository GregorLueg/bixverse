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

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

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

## Examples

``` r
# split the modules into positively and negatively correlated genes
mat <- t(synthetic_signal_matrix()$mat)
obj <- BulkCoExp(mat, data.table::data.table(sample_id = rownames(mat)))
obj <- preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
obj <- cor_module_processing(obj, cor_method = "spearman", .verbose = FALSE)
obj <- cor_module_coremo_clustering(obj, .verbose = FALSE)
obj <- cor_module_coremo_cor_sign(obj, .verbose = FALSE)
table(obj@outputs$final_modules$sign)
#> 
#> neg pos 
#>   2 218 
```
