# Run fgsea simple method for gene ontology with elimination method

Run fgsea simple method for gene ontology with elimination method

## Usage

``` r
rs_geom_elim_fgsea_simple(
  stats,
  levels,
  go_obj,
  gsea_params,
  elim_threshold,
  iters,
  seed
)
```

## Arguments

- stats:

  Named numerical vector. Needs to be sorted. The gene level statistics.

- levels:

  A character vector representing the levels to iterate through. The
  order will be the one the iterations are happening in.

- go_obj:

  The gene_ontology_data S7 class. See
  [`gene_ontology_data()`](gene_ontology_data.md).

- gsea_params:

  List. The GSEA parameters, see [`params_gsea()`](params_gsea.md)
  wrapper function. This function generates a list containing:

  - min_size - Integer. Minimum size for the gene sets.

  - max_size - Integer. Maximum size for the gene sets.

  - gsea_param - Float. The GSEA parameter. Defaults to `1.0`.

  - sample_size - Integer. Number of samples to iterate through for the
    multi-level implementation of fgsea.

  - eps - Float. Boundary for calculating the p-value. Used for the
    multi- level implementation of fgsea.

- elim_threshold:

  p-value below which the elimination procedure shall be applied to the
  ancestors.

- iters:

  Integer. Number of random permutations for the fgsea simple method to
  use

- seed:

  Integer. For reproducibility purposes.

## Value

List with the following elements

- go_ids The name of the tested gene ontology identifer.

- es The enrichment scores for the pathway

- nes The normalised enrichment scores for the pathway

- size The pathway sizes (after elimination!).

- pvals The p-values for this pathway based on permutation testing

- n_more_extreme Number of times the enrichment score was bigger or
  smaller than the permutation (pending sign).

- le_zero Number of times the permutation was less than zero.

- ge_zero Number of times the permutation was greater than zero.

- leading_edge A list of the index positions of the leading edge genes
  for this given GO term.
