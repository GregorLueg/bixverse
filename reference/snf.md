# Similarity network fusion (deprecated)

**\[deprecated\]**

This constructor has been renamed to
[`SimilarityNetworkFusion()`](https://gregorlueg.github.io/bixverse/reference/SimilarityNetworkFusion.md).

## Usage

``` r
snf(data = NULL, data_name = NULL, snf_params = params_snf())
```

## Arguments

- data:

  Optional data to transform into adjacency data. Can be a data.table
  (categorical/mixed) or a matrix (continuous).

- data_name:

  Optional string. Name of the data modality.

- snf_params:

  List. The SNF parameters, see
  [`params_snf()`](https://gregorlueg.github.io/bixverse/reference/params_snf.md).

## Value

Returns a
[`SimilarityNetworkFusion()`](https://gregorlueg.github.io/bixverse/reference/SimilarityNetworkFusion.md)
object.
