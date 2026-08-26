# Calculate VISION pathway scores in Rust (for meta cells)

**\[experimental\]** The function will take in a list of gene sets that
contains lists of `"pos"` and `"neg"` gene indices (0-indexed). You
don't have to provide the `"neg"`, but it can be useful to classify the
delta of two stats (EMT, Th1; Th2) etc. This version works on MetaCell
counts which are stored in memory directly.

## Usage

``` r
rs_mc_vision(sparse_data, gs_list, verbose)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `cs_type`. Shape is (metacells, genes) and the data need to
  be the **normalised** counts.

- gs_list:

  Nested list. Each sublist contains the (0-indexed!) positive and
  negative gene indices of that specific gene set.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A matrix of meta cells x vision scores per gene set.

## References

DeTomaso, et al., Nat. Commun., 2019
