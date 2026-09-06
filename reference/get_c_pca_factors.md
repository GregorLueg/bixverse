# Get the contrastive PCA factors

Getter function for the feature factors of the contrastive PCA

## Usage

``` r
get_c_pca_factors(object)
```

## Arguments

- object:

  The underlying class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

The sample scores of the contrastive PCA run. If not found, returns a
warning and NULL.

## Examples

``` r
# sample scores of the contrastive components
cpca_data <- synthetic_c_pca_data()
target <- t(cpca_data$target)
background <- t(cpca_data$background)
meta <- data.table::data.table(sample_id = rownames(target))
obj <- BulkCoExp(target, meta)
obj <- preprocess_bulk_coexp(obj, .verbose = FALSE)
obj <- contrastive_pca_processing(obj, background, .verbose = FALSE)
obj <- contrastive_pca(obj, alpha = 2.5, no_pcs = 5L)
head(get_c_pca_factors(obj)[, 1:3])
#>               cPC_1      cPC_2       cPC_3
#> sample_1  0.9703806  0.2352614  0.92219735
#> sample_2 -3.0559242  1.2606351 -0.78755725
#> sample_3  1.3357543  0.2804866  1.42740465
#> sample_4  0.6604632  0.2926076 -1.17078103
#> sample_5  0.3330629  0.6744638 -0.13867245
#> sample_6  0.6910951 -0.6427294 -0.02157507
```
