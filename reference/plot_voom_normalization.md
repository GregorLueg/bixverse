# Helper plot function for Voom normalisation

Helper plot function for Voom normalisation

## Usage

``` r
plot_voom_normalization(norm_counts)
```

## Arguments

- norm_counts:

  Numeric matrix. The voom log2-CPM values, genes x samples.

## Value

ggplot object, i.e., voom normalisation plot.

## Examples

``` r
# mean variance trend after voom
syn <- synthetic_bulk_cor_matrix()
grp <- rep(c("case", "control"), each = 50)
norm_counts <- rs_cpm(syn$counts, lib_size = NULL, log = TRUE,
  prior_count = 0.5)
plot_voom_normalization(norm_counts)
#> `geom_smooth()` using formula = 'y ~ x'
```
