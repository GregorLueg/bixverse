# Get the final SNF matrix

Get the final SNF matrix

## Usage

``` r
get_snf_final_mat(object)
```

## Arguments

- object:

  The underlying class
  [`SimilarityNetworkFusion()`](https://gregorlueg.github.io/bixverse/reference/SimilarityNetworkFusion.md).

## Value

Returns the SNF adjacency/similarity matrix.

## Examples

``` r
# fused similarity matrix across two modalities
set.seed(42)
continuous <- matrix(rnorm(120), nrow = 12, ncol = 10)
rownames(continuous) <- sprintf("sample_%02i", 1:12)
colnames(continuous) <- sprintf("feature_%i", 1:10)
clinical <- data.table::data.table(
  sample_id = rownames(continuous),
  sex = factor(sample(c("M", "F"), 12, replace = TRUE)),
  stage = factor(sample(c("I", "II", "III"), 12, replace = TRUE))
)
object <- SimilarityNetworkFusion(
  data = continuous,
  data_name = "continuous",
  snf_params = params_snf(k = 3L)
)
object <- add_snf_data_modality(object, clinical, data_name = "clinical")
object <- run_snf(object)
dim(get_snf_final_mat(object))
#> [1] 12 12
```
