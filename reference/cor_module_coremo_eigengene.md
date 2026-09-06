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

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

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

## Examples

``` r
# eigengenes per module and the gene to eigengene correlations
mat <- t(synthetic_signal_matrix()$mat)
obj <- BulkCoExp(mat, data.table::data.table(sample_id = rownames(mat)))
obj <- preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
obj <- cor_module_processing(obj, cor_method = "spearman", .verbose = FALSE)
obj <- cor_module_coremo_clustering(obj, .verbose = FALSE)
obj <- cor_module_coremo_eigengene(obj, .verbose = FALSE)
head(get_modules(get_results(obj)))
#> Key: <gene>
#>       gene module_id eigengene_cor
#>     <char>    <char>         <num>
#> 1:  gene10         1     0.7968888
#> 2: gene100         1     0.7727678
#> 3: gene101         2     0.5837753
#> 4: gene102         2     0.7664666
#> 5: gene103         2     0.7334704
#> 6: gene107         2     0.5509239
```
