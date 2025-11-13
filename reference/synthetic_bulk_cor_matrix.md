# Generates synthetic bulk RNAseq data

Function generates synthetic bulkRNAseq data with heteroskedasticity
(lowly expressed genes show higher variance) and can optionally add gene
~ gene correlations for testing purposes of module detection methods.

## Usage

``` r
synthetic_bulk_cor_matrix(
  num_samples = 100L,
  num_genes = 1000L,
  add_modules = TRUE,
  module_sizes = c(100L, 100L, 100L),
  seed = 123L
)
```

## Arguments

- num_samples:

  Integer. Number of samples.

- num_genes:

  Integer. Number of genes.

- add_modules:

  Boolean. Shall gene modules with correlation structures be generated.

- module_sizes:

  Integer vector. Sizes of the different correlation modules. The sum
  needs be smaller than `num_genes`.

- seed:

  Integer. Seed for reproducibility purposes.

## Value

A `synthetic_bulk_data` class containing:

- counts - The count matrix

- sparse_counts - A slot for sparse counts that can be added later, see
  [`simulate_dropouts()`](simulate_dropouts.md).

- module_data - The module membership of the genes
