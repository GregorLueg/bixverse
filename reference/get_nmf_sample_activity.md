# Get the NMF sample activity

Getter function to extract the sample activity matrix (samples x k) from
a bulk NMF fit. If no NMF fit is present, returns `NULL` with a warning.

## Usage

``` r
get_nmf_sample_activity(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

A samples x k numeric matrix (if found) or `NULL`.

## Examples

``` r
# sample activity of a four-factor fit
syn <- generate_gene_module_data(n_samples = 24L, n_genes = 60L)
# NMF needs a non-negative matrix
mat <- syn$data - min(syn$data)
obj <- BulkCoExp(mat, syn$meta_data)
obj <- preprocess_bulk_coexp(
  obj, hvg = NULL, scaling = FALSE, .verbose = FALSE
)
obj <- nmf_bulk(obj, k = 4L, .verbose = FALSE)
dim(get_nmf_sample_activity(obj))
#> [1] 24  4
```
