# Wrapper function to generate BulkDge object from h5ad

This is a helper function that can be used to create a `BulkDge` object
(see
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md))
directly from h5ad objects.

## Usage

``` r
bulk_dge_from_h5ad(h5_path, .verbose = TRUE)
```

## Arguments

- h5_path:

  String. Path to the h5ad object.

- .verbose:

  Controls verbosity of the function

## Value

`BulkDge` object.

## Examples

``` r
# round trip a synthetic count matrix through a temporary h5ad
syn <- synthetic_bulk_cor_matrix()
h5_path <- tempfile(fileext = ".h5ad")
write_h5ad_sc_dense(
  f_path = h5_path,
  counts = t(syn$counts),
  obs = data.table::data.table(sample_id = colnames(syn$counts)),
  var = data.table::data.table(var_id = rownames(syn$counts)),
  .verbose = FALSE
)
object <- bulk_dge_from_h5ad(h5_path, .verbose = FALSE)
object
#> Bulk differential gene expression class (BulkDge).
#>  Raw counts: 1000 genes x 100 samples.
#>  Meta-data rows: 100.
#>  Variable info provided: TRUE.
#>  Applied steps:
#>   qc_bulk_dge(): FALSE.
#>   normalise_bulk_dge(): FALSE.
#>   batch_correction_bulk_dge(): FALSE.
#>   calculate_pca_bulk_dge(): FALSE.
#>   calculate_dge_limma(): FALSE.
#>   calculate_dge_hedges(): FALSE.
#>   TPM normalisation: FALSE.
#>   FPKM normalisation: FALSE.

unlink(h5_path)
```
