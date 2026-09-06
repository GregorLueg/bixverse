# Run gene ontology enrichment with elimination method over a list.

This method takes the GeneOntologyElim and a list of target gene sets
and performs an GSE enrichment leveraging ontological information. It
starts at the lowest levels of the ontology and tests if there is
significant enrichment for any GO terms. If the threshold of the p-value
is below the elimination threshold, the genes from this term will be
removed from all its ancestors. The function then proceeds to the next
level of the ontology and repeats the process. The class will leverage
Rust threading to parallelise the process. The gene universe will be
automatically set to every gene represented in the ontology.

## Usage

``` r
gse_go_elim_method_list(
  object,
  target_gene_list,
  minimum_overlap = 3L,
  fdr_threshold = 0.05,
  elim_threshold = 0.05,
  min_genes = NULL
)
```

## Arguments

- object:

  The underlying class, see
  [`GeneOntologyElim()`](https://gregorlueg.github.io/bixverse/reference/GeneOntologyElim.md).

- target_gene_list:

  List. The target genes list you wish to apply the gene set enrichment
  analysis over.

- minimum_overlap:

  Integer. Threshold for the minimal overlap.

- fdr_threshold:

  Float. Threshold for maximum fdr to include in the output.

- elim_threshold:

  Float. Threshold from which p-value onwards the elimination on the
  ancestors shall be conducted.

- min_genes:

  Integer. Minimum number of genes that have to be included in the gene
  ontology term. If NULL, it will default to the number of minimum genes
  stored in `GeneOntologyElim`.

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
target_genes <- unique(unlist(S7::prop(go_obj, "go_to_genes")[1:5]))
res <- gse_go_elim_method_list(
  go_obj,
  target_gene_list = list(
    set_a = target_genes,
    set_b = rev(target_genes)[1:50]
  )
)
head(res, 3)
#>    target_set_name                                 go_name      go_id
#>             <char>                                  <char>     <char>
#> 1:           set_a                           axon guidance GO:0007411
#> 2:           set_a spliceosomal tri-snRNP complex assembly GO:0000244
#> 3:           set_b                           axon guidance GO:0007411
#>    odds_ratios         pvals           fdr  hits gene_set_lengths
#>          <num>         <num>         <num> <num>            <num>
#> 1:         Inf  0.000000e+00  0.000000e+00   192              192
#> 2:         Inf 2.939673e-140 2.850013e-137    76               76
#> 3:         Inf 1.104880e-105 2.155621e-102    50              192
# }
```
