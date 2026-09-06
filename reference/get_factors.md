# Get factor matrices from a BulkModuleResult

Returns one factor matrix by key, or the whole named list of factor
matrices if `which` is `NULL`. Keys are method-specific; see the
`factors` argument of
[`new_bulk_module_result()`](https://gregorlueg.github.io/bixverse/reference/new_bulk_module_result.md)
for the shared vocabulary.

## Usage

``` r
get_factors(object, which = NULL)
```

## Arguments

- object:

  A `BulkModuleResult`.

- which:

  String or `NULL`. Name of the factor matrix to return. If `NULL`,
  returns the whole list.

## Value

A matrix, the named list of matrices, or `NULL` (with warning) if
`which` is not among the stored factor keys.

## Examples

``` r
# gene loadings and sample activities of an NMF fit
syn <- synthetic_bulk_cor_matrix()
mat <- log1p(t(syn$counts))
meta <- data.table::data.table(sample_id = rownames(mat))
object <- BulkCoExp(raw_data = mat, meta_data = meta)
object <- preprocess_bulk_coexp(object, hvg = 500L, .verbose = FALSE)
object <- nmf_bulk(object, k = 3L, .verbose = FALSE)
res <- S7::prop(object, "final_results")
names(get_factors(res))
#> [1] "gene_loadings"   "sample_activity"
dim(get_factors(res, "gene_loadings"))
#> [1] 500   3
```
