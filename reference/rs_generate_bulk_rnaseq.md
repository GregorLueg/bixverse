# Generation of bulkRNAseq-like data with optional correlation structure

Function generates synthetic bulkRNAseq data with heteroskedasticity
(lowly expressed genes show higher variance) and can optionally add
correlation structures for testing purposes.

## Usage

``` r
rs_generate_bulk_rnaseq(
  num_samples,
  num_genes,
  seed,
  add_modules,
  module_sizes
)
```

## Arguments

- num_samples:

  Integer. Number of samples to simulate.

- num_genes:

  Integer. Number of genes to simulate.

- seed:

  Integer. Seed for reproducibility.

- add_modules:

  Boolean. Shall correlation structures be added to the data.

- module_sizes:

  `NULL` or vector of sizes of the gene modules. When `NULL` defaults to
  `c(300, 250, 200, 300, 500)`. Warning! The sum of this vector must be
  ≤ num_genes!

## Value

List with the following elements

- counts The matrix of simulated counts.

- module_membership Vector defining the module membership.
