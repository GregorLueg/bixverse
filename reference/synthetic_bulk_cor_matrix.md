# Generates synthetic bulk RNAseq data

Function generates synthetic bulkRNAseq data with heteroskedasticity
(lowly expressed genes show higher variance) and optional co-expression
modules for testing module detection methods. The ground truth comes
back with the counts, so a recovered module or eigengene can be scored
against what was actually simulated. Which generator to pick and what
each parameter does is documented in
[`params_synthetic_bulk_rnaseq()`](https://gregorlueg.github.io/bixverse/reference/params_synthetic_bulk_rnaseq.md).

## Usage

``` r
synthetic_bulk_cor_matrix(synthetic_params = params_synthetic_bulk_rnaseq())
```

## Arguments

- synthetic_params:

  List. The synthetic data parameters, see
  [`params_synthetic_bulk_rnaseq()`](https://gregorlueg.github.io/bixverse/reference/params_synthetic_bulk_rnaseq.md).

## Value

A `synthetic_bulk_data` class containing:

- counts - The count matrix. Rows are genes, columns are samples.

- sparse_counts - A slot for sparse counts that can be added later, see
  [`simulate_dropouts()`](https://gregorlueg.github.io/bixverse/reference/simulate_dropouts.md).

- module_data - Per-gene ground truth: the gene identifier, its
  `membership` (`0` for background), its `loading` on the module factor
  and whether it is a hub gene.

- module_factors - The latent factors the modules were built on. Rows
  are modules, columns are samples. This is what a module eigengene or
  an ICA/NMF component is trying to recover.

The parameters used are stored on the `synthetic_params` attribute.

## Examples

``` r
# three co-expression modules on a 1000 x 100 count matrix
syn <- synthetic_bulk_cor_matrix(params_synthetic_bulk_rnaseq())
dim(syn$counts)
#> [1] 1000  100
head(syn$module_data)
#>      gene membership   loading is_hub
#>    <char>      <int>     <num> <lgcl>
#> 1: gene_1          1 0.7308258  FALSE
#> 2: gene_2          1 5.5270005   TRUE
#> 3: gene_3          1 2.6595206   TRUE
#> 4: gene_4          1 1.0970747  FALSE
#> 5: gene_5          1 1.3758719  FALSE
#> 6: gene_6          1 0.6133037  FALSE
dim(syn$module_factors)
#> [1]   3 100
```
