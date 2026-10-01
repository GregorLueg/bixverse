# Build a Bonsai tree from metacell counts

**\[experimental\]** As
[`rs_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_bonsai.md),
with the metacells' aggregated raw counts from memory in place of the
binary files. Each metacell is a leaf, and its total counts over all
genes are its library size.

## Usage

``` r
rs_mc_bonsai(sparse_data, gene_indices, bonsai_params, verbose)
```

## Arguments

- sparse_data:

  List. The raw metacell counts, see
  [`mc_counts_to_list()`](https://gregorlueg.github.io/bixverse/reference/mc_counts_to_list.md)
  with `assay = "raw"`.

- gene_indices:

  Integer. The candidate genes. (0-indexed!)

- bonsai_params:

  List. Parameter list, see
  [`params_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/params_sc_bonsai.md).

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

The same list as
[`rs_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_bonsai.md),
with the metacells as the leaves in their row order.

## References

de Groot, et al., Nat Biotechnol, 2026; Breda, et al., Nat Biotechnol,
2021.
