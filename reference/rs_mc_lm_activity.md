# ULM / MLM activity scores (for meta cells)

**\[experimental\]** decoupler-style univariate (ULM) or multivariate
(MLM) linear model activity scores per meta cell on the normalised
counts, over all genes. This version works on MetaCell counts which are
stored in memory directly.

## Usage

``` r
rs_mc_lm_activity(sparse_data, gs_list, lm_params, verbose)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `cs_type`. Shape is (metacells, genes) and the data need to
  be the **normalised** counts.

- gs_list:

  List. One element per set, each a list with `idx` (integer, 0-indexed
  gene positions) and optionally `weight` (numeric, aligned with `idx`;
  absent means weight 1).

- lm_params:

  List. With `method` (`"ulm"` or `"mlm"`), `stat` (`"t"` or `"coef"`)
  and `p_values` (boolean).

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following elements

- es - Numerical matrix of meta cells x sets with the scores.

- p_values - Numerical matrix of meta cells x sets with the two-sided
  p-values, or `NULL` if not requested.

## References

Badia-i-Mompel, et al., Bioinform. Adv., 2022.
