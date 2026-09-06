# Prepare correlation-based module detection

This function will calculate the correlation coefficients between the
genes, using the highly variable genes (if available, otherwise the
function will use the raw data). The data will be stored in a
memory-efficient format in the properties of the class.

## Usage

``` r
cor_module_processing(
  object,
  cor_method = c("pearson", "spearman"),
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).
  Ideally, you should run
  [`preprocess_bulk_coexp()`](https://gregorlueg.github.io/bixverse/reference/preprocess_bulk_coexp.md)
  before applying this function.

- cor_method:

  String. Option of `c("pearson", "spearman")`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added data to the properties for subsequent usage.

## Examples

``` r
# spearman correlations over the 300 most variable genes
mat <- t(synthetic_signal_matrix()$mat)
obj <- BulkCoExp(mat, data.table::data.table(sample_id = rownames(mat)))
obj <- preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
obj <- cor_module_processing(obj, cor_method = "spearman", .verbose = FALSE)
dim(obj@processed_data$correlation_res$get_sym_matrix(
  .verbose = FALSE
))
#> [1] 300 300
```
