# Helper function to rapidly retrieve the indices of the gene set members

Helper function to rapidly retrieve the indices of the gene set members

## Usage

``` r
rs_get_gs_indices(gene_universe, pathway_list)
```

## Arguments

- gene_universe:

  Character Vector. The genes represented in the gene universe.

- pathway_list:

  List. A named list with each element containing the genes for this
  pathway.

## Value

Returns a list with the index positions of the gene set genes in the
gene universe. Importantly, these are indexed to R's 1-indexing!
