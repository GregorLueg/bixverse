# Gene set enrichment (GSE) based on a hypergeometric test over a list.

Takes a set of list of target genes, a list of gene sets and calculates
a p-value (hypergeometric test) and odds ratio (OR) against all the gene
sets. Also applies a multiple hypothesis correction (BH) to the
p-values.

## Usage

``` r
gse_hypergeometric_list(
  target_genes_list,
  gene_set_list,
  gene_universe = NULL,
  threshold = 0.05,
  minimum_overlap = 3L,
  .verbose = FALSE
)
```

## Arguments

- target_genes_list:

  Named list of character vectors. Names should represent the
  identifiers of the target genes and the elements the genes.

- gene_set_list:

  Named list of character vectors. Names should represent the gene sets,
  pathways, and the elements the genes within the respective gene set.

- gene_universe:

  Optional character vector. If you would like to specify specifically
  the gene universe. If set to NULL, the function will default to all
  represented genes in the `gene_set_list`.

- threshold:

  Float between 0 and 1 to filter on the fdr. Default: 0.05. If NULL
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
# two target sets tested against the same gene sets in one call
gene_universe <- sprintf("gene_%03i", 1:200)
gene_sets <- list(
  set_a = gene_universe[1:20],
  set_b = gene_universe[100:130]
)
targets <- list(
  hit_a = gene_universe[c(1:12, 150:158)],
  hit_b = gene_universe[c(100:112, 5:8)]
)
gse_hypergeometric_list(targets, gene_sets, gene_universe, threshold = 1)
#>    target_set_name odds_ratios        pvals          fdr  hits gene_set_lengths
#>             <char>       <num>        <num>        <num> <num>            <num>
#> 1:           hit_b   29.791667 3.811374e-09 7.622748e-09    13               31
#> 2:           hit_a   28.500000 4.196950e-09 8.393901e-09    12               20
#> 3:           hit_b    3.211538 7.369238e-02 7.369238e-02     4               20
#>    gene_set_name target_set_lengths
#>           <char>              <int>
#> 1:         set_b                 17
#> 2:         set_a                 21
#> 3:         set_a                 17
```
