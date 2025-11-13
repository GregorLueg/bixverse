# Bixverse implementation of the simple fgsea algorithm

Rust-based version of the fgsea simple algorithm. This one is
permutation- based and similar to the traditional implementation, but
leverages some clear tricks to be way faster, see Korotkevich, et al.

## Usage

``` r
calc_fgsea_simple(
  stats,
  pathways,
  nperm = 2000L,
  gsea_params = params_gsea(),
  seed = 123L
)
```

## Arguments

- stats:

  Named numeric vector. The gene level statistic.

- pathways:

  List. A named list with each element containing the genes for this
  pathway.

- nperm:

  Integer. Number of permutation tests. Defaults to `2000L`.

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

- seed:

  Random seed for reproducibility.

## Value

A data.table with the results from the GSEA with the following columns:

- es - Float. The enrichment score for this pathway.

- nes - Float. The normalised enrichment score for this pathway.

- pvals - Float. The p-value for this pathway.

- n_more_extreme - Integer. Number of permutation that had more extreme
  enrichment scores than the actual.

- size - Integer. The size of the pathway.

- pathway_name - Character. The name of the pathway.

- leading_edge - List of character vectors with the leading edge genes.

- fdr - Float. The adjusted pval.

## References

Korotkevich, et al., bioRxiv
