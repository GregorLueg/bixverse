# Run gene ontology enrichment with elimination method.

This method takes the GeneOntologyElim and a target gene set and
performs an GSE enrichment leveraging ontological information. It starts
at the lowest levels of the ontology and tests if there is significant
enrichment for any GO terms. If the threshold of the p-value is below
the elimination threshold, the genes from this term will be removed from
all its ancestors. The function then proceeds to the next level of the
ontology and repeats the process. The gene universe will be
automatically set to every gene represented in the ontology.

## Usage

``` r
gse_go_elim_method(
  object,
  target_genes,
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

- target_genes:

  String. The target genes you wish to apply the GSEA over.

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
head(gse_go_elim_method(go_obj, target_genes = target_genes), 3)
#>         go_id                                  go_name odds_ratios
#>        <char>                                   <char>       <num>
#> 1: GO:0007411                            axon guidance         Inf
#> 2: GO:0000244  spliceosomal tri-snRNP complex assembly         Inf
#> 3: GO:0000353 formation of quadruple SL/U4/U5/U6 snRNP         Inf
#>            pvals           fdr  hits gene_set_lengths
#>            <num>         <num> <num>            <num>
#> 1:  0.000000e+00  0.000000e+00   192              192
#> 2: 2.939673e-140 2.850013e-137    76               76
#> 3: 1.335653e-104 8.632771e-102    57               57
# }
```
