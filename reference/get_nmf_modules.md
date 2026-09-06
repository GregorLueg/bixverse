# Get the NMF module membership data.table

Getter function to extract the gene-to-module data.table from a bulk NMF
fit. Each row is one gene assigned to its top-loading module.

## Usage

``` r
get_nmf_modules(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

A data.table with `gene`, `module_id`, `loading`, `sign` columns (if
found) or `NULL`.

## Examples

``` r
# gene to module assignments of a four-factor fit
syn <- generate_gene_module_data(n_samples = 24L, n_genes = 60L)
# NMF needs a non-negative matrix
mat <- syn$data - min(syn$data)
obj <- BulkCoExp(mat, syn$meta_data)
obj <- preprocess_bulk_coexp(
  obj, hvg = NULL, scaling = FALSE, .verbose = FALSE
)
obj <- nmf_bulk(obj, k = 4L, .verbose = FALSE)
head(get_nmf_modules(obj))
#>          gene module_id  loading   sign        z
#>        <char>    <char>    <num> <char>    <num>
#> 1: feature_44   comp_01 3.675330    pos 38.28581
#> 2: feature_40   comp_01 3.643421    pos 37.94433
#> 3: feature_42   comp_01 3.633641    pos 37.83967
#> 4: feature_31   comp_01 3.630990    pos 37.81130
#> 5: feature_34   comp_01 3.602253    pos 37.50377
#> 6: feature_45   comp_01 3.590235    pos 37.37517
```
