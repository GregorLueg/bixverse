# Run the SNF algorithm

This function will run the SNF algorithm on top of the adjacency
matrices found in the object. You can also optionally specify which
adjacency matrices to use via the `to_include` parameter.

## Usage

``` r
run_snf(object, to_include = NULL, params = NULL)
```

## Arguments

- object:

  The underlying class, see
  [`SimilarityNetworkFusion()`](https://gregorlueg.github.io/bixverse/reference/SimilarityNetworkFusion.md).

- to_include:

  Optional string, if you wish to only use a subset of the generated
  adjacency matrices. If `NULL` all matrices will be used for the fusion
  process.

- params:

  Optional List. If you wish to overwite the already set up parameters
  for SNF, see
  [`params_snf()`](https://gregorlueg.github.io/bixverse/reference/params_snf.md).
  If `NULL`, the settings from within the object will be used. If not
  NULL, the new parameters will be used for this modality specifically
  and only for this modality!

## Value

The class with added adjacency matrix based on the SNF algorithm.

## Examples

``` r
# fuse a continuous and a categorical modality
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
