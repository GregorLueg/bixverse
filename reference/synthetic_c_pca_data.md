# Generates synthetic data for contrastive PCA exploration.

Generates three elements: target matrix, background matrix and labels
for the target matrix.

## Usage

``` r
synthetic_c_pca_data(seed = 10101L)
```

## Arguments

- seed:

  Integer. Initial random seed for generation of the synthetic data.
  Default: 10101L.

## Value

A `cpca_synthetic_data` class with the following elements:

- target - The target matrix.

- background - The background matrix.

- target_labels - The target labels

## Examples

``` r
# target and background matrices for a contrastive PCA run
cpca_data <- synthetic_c_pca_data()
dim(cpca_data$target)
#> [1]  30 400
dim(cpca_data$background)
#> [1]  30 400
table(cpca_data$target_labels)
#> 
#>   0   1   2   3 
#> 100 100 100 100 
```
