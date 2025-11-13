# Run gene ontology enrichment with elimination method over a list.

This method takes the gene_ontology_data and a list of target gene sets
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
  [`gene_ontology_data()`](gene_ontology_data.md).

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
  stored in `gene_ontology_data`.

## Value

data.table with enrichment results.
