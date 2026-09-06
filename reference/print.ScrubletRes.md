# Print a ScrubletRes object

Print a ScrubletRes object

## Usage

``` r
# S3 method for class 'ScrubletRes'
print(x, ...)
```

## Arguments

- x:

  A `ScrubletRes` object.

- ...:

  Ignored.

## Value

Invisible `x`.

## Examples

``` r
# the doublet rates and the threshold that produced them
sc <- demo_single_cells(prepped = FALSE)
res <- scrublet_sc(
  sc,
  scrublet_params = params_scrublet(
    pca = list(no_pcs = 10L),
    hvg = list(min_gene_var_pctl = 0.0),
    n_bins = 20L
  ),
  .verbose = FALSE
)
print(res)
#> ScrubletRes: 500 cells, 36 doublets (7.2%)
#>   Threshold:              0.1215
#>   Detected doublet rate:  7.2%
#>   Detectable fraction:    97.7%
#>   Overall doublet rate:   7.4%
#>   Simulated doublets:     750

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
