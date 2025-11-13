# Run a hypergeometric test over a list of target genes

Given a list of target gene sets, this function will test for each of
the individual target genes the hypergeoemetric enrichment against the
specified gene sets. WARNING! Incorrect use can cause kernel crashes.
Wrapper around the Rust functions with type checks are provided in the
package.

## Usage

``` r
rs_hypergeom_test_list(
  target_genes_list,
  gene_sets,
  gene_universe,
  min_overlap,
  fdr_threshold
)
```

## Arguments

- target_genes_list:

  A character vector representing the target gene set.

- gene_sets:

  A list of strings that represent the gene sets to test against.

- gene_universe:

  A character vector representing the gene universe from which the
  target genes and gene sets are sampled from.

- min_overlap:

  Optional integer. Shall a filter be applied on the minimum of
  overlappign genes.

- fdr_threshold:

  Optional float. Shall a filter be applied for the maximum tolerated
  FDR.

## Value

A list containing:

- pvals - The p-values from the hypergeometric test.

- fdr - The FDRs for each target gene calculated across all gene sets.

- odds ratios - The calculated odds ratios

- hits - The size of the overlap between the target gene set and
  individual gene sets.

- gene_set_lengths - The length of the gene sets.

- to_keep - Indices of the tests that passed.

- tests_passed - How many tests passed the filter criteria for that
  target set.
