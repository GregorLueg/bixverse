# Calculate AUCell in Rust (for meta cells)

**\[experimental\]** The function will take in a list of gene set
indices (0-indexed!) and calculate an AUCell type statistic. Three
options here: the recovery-curve AUC of Aibar, et al. (the actual AUCell
statistic), an AUC derived from the Mann-Whitney statistic, or average
precision. This version works on MetaCell counts which are stored in
memory directly.

## Usage

``` r
rs_mc_aucell(sparse_data, gs_list, aucell_params, verbose)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `format`.

- gs_list:

  List. List with the gene set indices (0-indexed!) of the genes of
  interest.

- aucell_params:

  List. The AUCell parameters, see
  [`params_sc_aucell()`](https://gregorlueg.github.io/bixverse/reference/params_sc_aucell.md).

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A matrix of cells x gene sets with the values representing the AUC.
