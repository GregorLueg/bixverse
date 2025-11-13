# Generate reciprocal best hits based on set similarities

This function takes a nested list that contains gene modules/ sets
derived from various methods and generate identifies reciprocal best
hits between gene modules/sets across the different origins. WARNING!
Incorrect use can cause kernel crashes. Wrapper around the Rust
functions with type checks are provided in the package.

## Usage

``` r
rs_rbh_sets(module_list, overlap_coefficient, min_similarity)
```

## Arguments

- module_list:

  A nested named list. The outer list should contain the origin of the
  gene modules, the inner list the names of the gene modules and the
  respective genes in them.

- overlap_coefficient:

  Shall the overlap coefficient instead of the Jaccard similarity be
  used.

- min_similarity:

  Minimum similarity that should exist between any two given gene
  modules to actually calculate RBH pairs.

## Value

A list containing:

- origin - The name of the origin of the gene modules.

- target - The name of the target of the gene modules.

- comparisons - Integer vector indicating how many RBH hits were
  identified in this comparison

- origin_modules - Names of the gene modules from the origin.

- target_modules - Names of the gene modules from the target.

- similarity - The similarities between the two respective gene modules.
