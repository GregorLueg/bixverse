# Get an individual affinity matrix

Get an individual affinity matrix

## Usage

``` r
get_snf_adjcacency_mat(object, name)
```

## Arguments

- object:

  The underlying class
  [`SimilarityNetworkFusion()`](https://gregorlueg.github.io/bixverse/reference/SimilarityNetworkFusion.md).

- name:

  String. The name of the individual data modality affinity matrix to
  return.

## Value

Returns adjcacency matrix if found.

## Examples

``` r
# affinity matrix of a single modality
set.seed(42)
continuous <- matrix(rnorm(120), nrow = 12, ncol = 10)
rownames(continuous) <- sprintf("sample_%02i", 1:12)
colnames(continuous) <- sprintf("feature_%i", 1:10)
object <- SimilarityNetworkFusion(
  data = continuous,
  data_name = "continuous",
  snf_params = params_snf(k = 3L)
)
dim(get_snf_adjcacency_mat(object, name = "continuous"))
#> [1] 12 12
```
