# Assesses CoReMo-based gene module stability

The function assesses the stability of the CoReMo modules, leveraging a
leave-one-out sample method. In each iteration, one of the samples is
left out and the clustering is redone. Subsequently, the Jaccard
similarity (as a surrogate for membership stability) is calculated for
each feature across the different resamplings and added to the
`final_modules` data.table.

## Usage

``` r
cor_module_coremo_stability(object, chunk_size = 15L, .verbose = TRUE)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

- chunk_size:

  Integer. Chunk size in which to process the data. Defaults to `15L`,
  i.e., 15 samples are being processed in one go. You can use bigger
  values here, but be aware of memory pressure.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added data to the properties for subsequent usage.

## References

Srivastava, et al., Nat. Commun., 2018; Francois, Romagnolo, et al.,
Nat. Commun., 2024.

## Examples

``` r
# leave-one-out stability of the CoReMo modules
mat <- t(synthetic_signal_matrix()$mat)
obj <- BulkCoExp(mat, data.table::data.table(sample_id = rownames(mat)))
obj <- preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
obj <- cor_module_processing(obj, cor_method = "spearman", .verbose = FALSE)
obj <- cor_module_coremo_clustering(obj, .verbose = FALSE)
obj <- cor_module_coremo_stability(obj, .verbose = FALSE)
summary(obj@outputs$final_modules$stability)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.3035  0.8156  0.8771  0.8065  0.9041  0.9062 
```
