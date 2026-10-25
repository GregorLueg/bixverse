# Wrapper function for the MetaCells2 divide and conquer parameters

Wrapper function for the MetaCells2 divide and conquer parameters

## Usage

``` r
params_sc_metacells2_dac(
  min_target_pile_size = 8000L,
  max_target_pile_size = 24000L,
  target_metacells_in_pile = 100L,
  target_pile_size = NULL,
  quick_and_dirty = FALSE,
  max_parallel_piles = 1L,
  groups_similarity_method = c("abs_pearson", "pearson"),
  rare_gene_modules = TRUE,
  rare_max_gene_cell_fraction = 0.001,
  rare_max_genes = 500L,
  rare_min_genes_of_modules = 4L,
  rare_min_cells_of_modules = 12L,
  rare_genes_similarity_method = c("repeated_pearson", "repeated_abs_pearson", "pearson",
    "abs_pearson"),
  rare_genes_cluster_method = c("ward", "complete", "average", "weighted")
)
```

## Arguments

- min_target_pile_size:

  Integer. Lower clamp on the derived pile size. Defaults to `8000L`.

- max_target_pile_size:

  Integer. Upper clamp on the derived pile size. Defaults to `24000L`.

- target_metacells_in_pile:

  Integer. Meta cells wanted per pile, used to derive the pile size.
  Defaults to `100L`.

- target_pile_size:

  Integer or `NULL`. Fixed pile size. `NULL` derives it from the data
  and clamps it to `[min_target_pile_size, max_target_pile_size]`.
  Defaults to `NULL`.

- quick_and_dirty:

  Boolean. Stop after the final phase instead of regrouping the
  outliers. Defaults to `FALSE`.

- max_parallel_piles:

  Integer. Piles computed concurrently. Each pile holds two dense
  `pile_size x pile_size` matrices, so this trades memory for speed.
  Defaults to `1L`.

- groups_similarity_method:

  String. Similarity used when grouping meta cells into piles. One of
  `c("abs_pearson", "pearson")`. Defaults to `"abs_pearson"`.

- rare_gene_modules:

  Boolean. Detect rare gene modules first and pile their cells
  separately. Defaults to `TRUE`.

- rare_max_gene_cell_fraction:

  Numeric. Largest fraction of cells a gene may be expressed in to be a
  rare gene candidate. Defaults to `0.001`.

- rare_max_genes:

  Integer. Most rare gene candidates kept. Defaults to `500L`.

- rare_min_genes_of_modules:

  Integer. Fewest genes in a rare gene module. Defaults to `4L`.

- rare_min_cells_of_modules:

  Integer. Fewest cells in a rare gene module. Defaults to `12L`.

- rare_genes_similarity_method:

  String. Gene-gene similarity of the rare candidates. One of
  `c("repeated_pearson", "repeated_abs_pearson", "pearson", "abs_pearson")`.
  Defaults to `"repeated_pearson"`.

- rare_genes_cluster_method:

  String. Linkage for the hierarchical clustering of the rare
  candidates. One of `c("ward", "complete", "average", "weighted")`.
  Defaults to `"ward"`.

## Value

A named list with the following elements:

- min_target_pile_size - Integer. Lower clamp on the derived pile size.
  Defaults to `8000L`.

- max_target_pile_size - Integer. Upper clamp on the derived pile size.
  Defaults to `24000L`.

- target_metacells_in_pile - Integer. Meta cells wanted per pile, used
  to derive the pile size. Defaults to `100L`.

- target_pile_size - Integer or `NULL`. Fixed pile size. `NULL` derives
  it from the data and clamps it to
  `[min_target_pile_size, max_target_pile_size]`. Defaults to `NULL`.

- quick_and_dirty - Boolean. Stop after the final phase instead of
  regrouping the outliers. Defaults to `FALSE`.

- max_parallel_piles - Integer. Piles computed concurrently. Each pile
  holds two dense `pile_size x pile_size` matrices, so this trades
  memory for speed. Defaults to `1L`.

- groups_similarity_method - String. Similarity used when grouping meta
  cells into piles. One of `c("abs_pearson", "pearson")`. Defaults to
  `"abs_pearson"`.

- rare_gene_modules - Boolean. Detect rare gene modules first and pile
  their cells separately. Defaults to `TRUE`.

- rare_max_gene_cell_fraction - Numeric. Largest fraction of cells a
  gene may be expressed in to be a rare gene candidate. Defaults to
  `0.001`.

- rare_max_genes - Integer. Most rare gene candidates kept. Defaults to
  `500L`.

- rare_min_genes_of_modules - Integer. Fewest genes in a rare gene
  module. Defaults to `4L`.

- rare_min_cells_of_modules - Integer. Fewest cells in a rare gene
  module. Defaults to `12L`.

- rare_genes_similarity_method - String. Gene-gene similarity of the
  rare candidates. One of
  `c("repeated_pearson", "repeated_abs_pearson", "pearson", "abs_pearson")`.
  Defaults to `"repeated_pearson"`.

- rare_genes_cluster_method - String. Linkage for the hierarchical
  clustering of the rare candidates. One of
  `c("ward", "complete", "average", "weighted")`. Defaults to `"ward"`.

## References

Ben-Kiki, et al., Genome Biol., 2022
