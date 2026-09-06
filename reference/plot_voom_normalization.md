# Helper plot function for Voom normalisation

Helper plot function for Voom normalisation

## Usage

``` r
plot_voom_normalization(voom_object)
```

## Arguments

- voom_object:

  `EList`. Voom object with normalised counts.

## Value

ggplot object, i.e., voom normalisation plot.

## Examples

``` r
# mean variance trend after voom
syn <- synthetic_bulk_cor_matrix()
grp <- rep(c("case", "control"), each = 50)
dge_list <- edgeR::normLibSizes(edgeR::DGEList(counts = syn$counts))
voom_obj <- limma::voom(dge_list, stats::model.matrix(~grp))
plot_voom_normalization(voom_obj)
#> `geom_smooth()` using formula = 'y ~ x'
```
