# Run GO enrichment with elimination with fgsea simple

This method takes the GeneOntologyElim and a vector of gene level
statistics to perform fgsea (simple) leveraging ontological information.
It starts at the lowest levels of the ontology and tests if there is
significant enrichment for any GO terms. If the threshold of the p-value
is below the elimination threshold, the genes from this term will be
removed from all its ancestors. The function then proceeds to the next
level of the ontology and repeats the process.

## Usage

``` r
fgsea_simple_go_elim(
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
res <- fgsea_simple_go_elim(go_obj, stats = stats, nperm = 1000L)
head(data.table::setorder(res, pvals), 3)
#>         go_id        es      nes  size       pvals
#>        <char>     <num>    <num> <num>       <num>
#> 1: GO:0007411 0.4597447 2.227786   192 0.001742160
#> 2: GO:0005685 0.7636759 3.059545    62 0.001795332
#> 3: GO:0000244 0.8595015 3.597023    76 0.001798561
#>                                                                                               leading_edge
#>                                                                                                     <list>
#> 1: ENSG00000274428,ENSG00000010810,ENSG00000278048,ENSG00000199568,ENSG00000287979,ENSG00000125753,...[72]
#> 2: ENSG00000011201,ENSG00000082516,ENSG00000177706,ENSG00000283527,ENSG00000206052,ENSG00000206625,...[41]
#> 3: ENSG00000196189,ENSG00000119335,ENSG00000164040,ENSG00000085552,ENSG00000169306,ENSG00000288093,...[67]
#>                                    go_name       fdr
#>                                     <char>     <num>
#> 1:                           axon guidance 0.2426716
#> 2:                                U1 snRNP 0.2426716
#> 3: spliceosomal tri-snRNP complex assembly 0.2426716
# }
```
