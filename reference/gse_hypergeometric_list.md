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
