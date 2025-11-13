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
