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

- meta_data - The sample metadata.

## Examples

``` r
# four partially overlapping modules over 24 samples
mods <- generate_gene_module_data(
  n_samples = 24L,
  n_genes = 60L,
  n_modules = 4L
)
dim(mods$data)
#> [1] 24 60
head(mods$meta_data)
#>    sample_id module_1_active module_2_active module_3_active module_4_active
#>       <char>          <lgcl>          <lgcl>          <lgcl>          <lgcl>
#> 1:  sample_1            TRUE           FALSE           FALSE           FALSE
#> 2:  sample_2            TRUE           FALSE           FALSE           FALSE
#> 3:  sample_3            TRUE           FALSE           FALSE           FALSE
#> 4:  sample_4            TRUE           FALSE           FALSE           FALSE
#> 5:  sample_5            TRUE            TRUE           FALSE           FALSE
#> 6:  sample_6            TRUE            TRUE           FALSE           FALSE
```
