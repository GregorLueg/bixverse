# Run hypergeometric enrichment over the gene ontology

**\[experimental\]** This function implements a Rust version of the gene
ontology enrichment with elimination: the starting point are the leafs
of the ontology and hypergeometric tests will first conducted there.
Should the hypergeometric test p-value be below a certain threshold, the
genes of that gene ontology term will be removed from all ancestors.

## Usage

``` r
rs_gse_geom_elim(
  target_genes,
  levels,
  go_obj,
  gene_universe_length,
  min_genes,
  elim_threshold,
  min_overlap,
  fdr_threshold
)
```

## Arguments

- target_genes:

  A character vector representing the target gene set.

- levels:

  A character vector representing the levels to iterate through. The
  order will be the one the iterations are happening in.

- go_obj:

  The `GeneOntologyElim` S7 class. See
  [`GeneOntologyElim()`](https://gregorlueg.github.io/bixverse/reference/GeneOntologyElim.md).

- gene_universe_length:

  The length of the gene universe.

- min_genes:

  number of minimum genes for the gene ontology term to be tested.

- elim_threshold:

  p-value below which the elimination procedure shall be applied to the
  ancestors.

- min_overlap:

  Optional minimum overlap threshold.

- fdr_threshold:

  Optional fdr threshold.

## Value

A list containing:

- go_ids - The gene ontology identifier.

- pvals - The calculated odds ratios.

- odds_ratios - The calculated odds ratios.

- overlap - The size of the overlap.

- gene_set_lengths - The length of the gene sets.
