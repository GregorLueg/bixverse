# TPM calculation

TPM calculation

## Usage

``` r
calculate_tpm(counts, gene_lengths)
```

## Arguments

- counts:

  Numeric matrix. Count matrix (gene x sample)

- gene_lengths:

  Named vector. Named vector with gene lengths.

## Value

TPM-normalised matrix.

## Examples

``` r
# TPM over synthetic counts with flat 2kb gene lengths
syn <- synthetic_bulk_cor_matrix()
gene_lengths <- stats::setNames(
  rep(2000, nrow(syn$counts)),
  rownames(syn$counts)
)
tpm <- calculate_tpm(syn$counts, gene_lengths)
colSums(tpm)[1:3]
#> sample_1 sample_2 sample_3 
#>    1e+06    1e+06    1e+06 
```
