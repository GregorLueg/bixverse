# Get the NMF gene loadings

Getter function to extract the gene loadings matrix (features x k) from
a bulk NMF fit stored in
[`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).
If no NMF fit is present, returns `NULL` with a warning.

## Usage

``` r
get_nmf_gene_loadings(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

A features x k numeric matrix (if found) or `NULL`.

## Examples

``` r
# gene loadings of a four-factor fit
syn <- generate_gene_module_data(n_samples = 24L, n_genes = 60L)
# NMF needs a non-negative matrix
mat <- syn$data - min(syn$data)
obj <- BulkCoExp(mat, syn$meta_data)
obj <- preprocess_bulk_coexp(
  obj, hvg = NULL, scaling = FALSE, .verbose = FALSE
)
obj <- nmf_bulk(obj, k = 4L, .verbose = FALSE)
dim(get_nmf_gene_loadings(obj))
#> [1] 60  4
```
