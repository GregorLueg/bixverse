# Run GO enrichment with elimination method over a continuous vectors

This method takes the GeneOntologyElim and a vector of gene level
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
  [`GeneOntologyElim()`](https://gregorlueg.github.io/bixverse/reference/GeneOntologyElim.md).

- stats:

  Named numeric vector. The gene level statistic.

- elim_threshold:

  Float. Threshold from which p-value onwards the elimination on the
  ancestors shall be conducted.

- nperm:

  Integer. Number of permutation tests. Defaults to `2000L`

- gsea_params:

  List. The GSEA parameters, see
  [`params_gsea()`](https://gregorlueg.github.io/bixverse/reference/params_gsea.md)
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

## Examples

``` r
# \donttest{
# human GO terms with at least 25 genes
go_obj <- GeneOntologyElim(
  get_go_data_human(.verbose = FALSE),
  min_genes = 25L
)
genes <- unique(unlist(S7::prop(go_obj, "go_to_genes")))
set.seed(1L)
stats <- stats::setNames(rnorm(length(genes)), genes)
stats[1:200] <- stats[1:200] + 2
res <- fgsea_go_elim(go_obj, stats = stats, nperm = 1000L)
head(res, 3)
#>         go_id        es      nes  size        pvals n_more_extreme
#>        <char>     <num>    <num> <num>        <num>          <num>
#> 1: GO:0000244 0.8595015 3.597023    76 6.153695e-26              0
#> 2: GO:0046540 0.7694975 3.312726    94 6.220826e-26              0
#> 3: GO:0000353 0.8578760 3.360251    57 6.254944e-26              0
#>                                                                                               leading_edge
#>                                                                                                     <list>
#> 1: ENSG00000199568,ENSG00000283527,ENSG00000275174,ENSG00000206625,ENSG00000288093,ENSG00000101161,...[67]
#> 2: ENSG00000199568,ENSG00000283527,ENSG00000275174,ENSG00000206625,ENSG00000288093,ENSG00000101161,...[62]
#> 3: ENSG00000199568,ENSG00000283527,ENSG00000275174,ENSG00000206625,ENSG00000288093,ENSG00000283372,...[50]
#>    log2err          fdr                                  go_name
#>      <num>        <num>                                   <char>
#> 1:      NA 3.080252e-23  spliceosomal tri-snRNP complex assembly
#> 2:      NA 3.080252e-23             U4/U6 x U5 tri-snRNP complex
#> 3:      NA 3.080252e-23 formation of quadruple SL/U4/U5/U6 snRNP
# }
```
