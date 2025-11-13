# Run GO enrichment with elimination method over a continuous vectors

This method takes the gene_ontology_data and a vector of gene level
statistics to perform fgsea (multi-level) leveraging ontological
information. It starts at the lowest levels of the ontology and tests if
there is significant enrichment for any GO terms. If the threshold of
the p-value is below the elimination threshold, the genes from this term
will be removed from all its ancestors. The function then proceeds to
the next level of the ontology and repeats the process. Subsequently, it
leverages the multi-level method to estimate lower p-values for
significant terms, see Korotkevich, et al.

## Usage

``` r
fgsea_go_elim(
  object,
  stats,
  elim_threshold = 0.05,
  nperm = 2000L,
  gsea_params = params_gsea(max_size = 2000L),
  seed = 42L
)
```

## Arguments

- object:

  The underlying class, see
  [`gene_ontology_data()`](gene_ontology_data.md).

- stats:

  Named numeric vector. The gene level statistic.

- elim_threshold:

  Float. Threshold from which p-value onwards the elimination on the
  ancestors shall be conducted.

- nperm:

  Integer. Number of permutation tests. Defaults to `2000L`

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

data.table with enrichment results.

## References

Korotkevich, et al., bioRxiv
