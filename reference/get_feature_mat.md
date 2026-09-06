# Get the feature matrix used for the classifier

Get the feature matrix used for the classifier

## Usage

``` r
get_feature_mat(x, ...)

# S3 method for class 'ScDblFinderRes'
get_feature_mat(x, ...)
```

## Arguments

- x:

  An object to get the feature matrix from. This will only include the
  values of the observed cells.

- ...:

  Additional parameters to forward to the method.

## Examples

``` r
# the features the scDblFinder classifier was trained on
sc <- demo_single_cells(prepped = FALSE)
res <- scdblfinder_sc(
  sc,
  scdblfinder_params = params_scdblfinder(
    pca = list(no_pcs = 10L),
    n_genes = 25L,
    cxds_genes = 25L
  ),
  return_features = TRUE,
  .verbose = FALSE
)
dim(get_feature_mat(res))
#> [1] 500  22

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
