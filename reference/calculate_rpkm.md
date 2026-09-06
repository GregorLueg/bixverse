# RPKM calculation

RPKM calculation

## Usage

``` r
calculate_rpkm(counts, gene_lengths)
```

## Arguments

- counts:

  Numeric matrix. Count matrix (gene x sample)

- gene_lengths:

  Named vector. Named vector with gene lengths.

## Value

RPKM-normalised matrix.

## Examples

``` r
# RPKM over synthetic counts with flat 2kb gene lengths
syn <- synthetic_bulk_cor_matrix()
gene_lengths <- stats::setNames(
  rep(2000, nrow(syn$counts)),
  rownames(syn$counts)
)
rpkm <- calculate_rpkm(syn$counts, gene_lengths)
rpkm[1:3, 1:3]
#>        sample_1  sample_2  sample_3
#> gene_1 363.1402  713.9413 339.18333
#> gene_2   0.0000  315.2468  54.70699
#> gene_3  22.0085 1019.9162 339.18333
```
