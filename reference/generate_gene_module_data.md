# Generates synthetic gene module data.

Generates an artifical matrix of active modules (for a maximum of 8) in
samples being active. Also allows for overlap between the modules.
Designed to test some of the algorithms.

## Usage

``` r
generate_gene_module_data(
  n_samples = 24L,
  n_genes = 60L,
  n_modules = 4L,
  overlap_fraction = 0.25,
  seed = 10101L
)
```

## Arguments

- n_samples:

  Integer. Number of samples.

- n_genes:

  Integer. Number of genes.

- n_modules:

  Integer. Number of active gene modules. To a maximum of 10.

- overlap_fraction:

  Float. How much the same modules can be active in the same samples.

- seed:

  Integer. Initial random seed for generation of the synthetic data.
  Default: 10101L.

## Value

A `synthetic_matrix_modules` class with the following items:

- data - The data matrix.

- metadata - The sample metadata.
