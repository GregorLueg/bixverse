# Gene set enrichment (GSE) based on a hypergeometric test.

Takes a set of target genes, a list of gene sets and calculates a
p-value (hypergeometric test) and odds ratio (OR) against all the gene
sets. Also applies a multiple hypothesis correction (BH) to the
p-values.

## Usage

``` r
gse_hypergeometric(
  target_genes,
  gene_set_list,
  gene_universe = NULL,
  threshold = 0.05,
  minimum_overlap = 3L,
  .verbose = FALSE
)
```

## Arguments

- target_genes:

  Character vector. GeneID(s) of the target genes.

- gene_set_list:

  Named list of character vectors. Names should represent the gene sets,
  pathways, and the elements the genes within the respective gene set.

- gene_universe:

  Optional character vector. If you would like to specify specifically
  the gene universe. If set to NULL, the function will default to all
  represented genes in the `gene_set_list`.

- threshold:

  Float between 0 and 1 to filter on the fdr. Default: 0.05. If 1
  everything is returned.

- minimum_overlap:

  Number of minimum overlap between the target genes and the respective
  gene set.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

data.table with enrichment results.

## Examples

``` r
# hypergeometric test of a target set against a small universe
gene_universe <- sprintf("gene_%03i", 1:200)
gene_sets <- list(
  set_a = gene_universe[1:20],
  set_b = gene_universe[15:40],
  set_c = gene_universe[100:130]
)
target <- gene_universe[c(1:12, 150:158)]
gse_hypergeometric(target, gene_sets, gene_universe, threshold = 1)
#>    gene_set_name odds_ratios       pvals          fdr  hits gene_set_lengths
#>           <char>       <num>       <num>        <num> <num>            <num>
#> 1:         set_a        28.5 4.19695e-09 1.259085e-08    12               20
#>    target_set_lengths
#>                 <int>
#> 1:                 21
```
