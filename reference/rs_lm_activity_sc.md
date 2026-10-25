# ULM / MLM activity scores per cell (decoupler-style)

**\[experimental\]** Fits a univariate (ULM) or multivariate (MLM)
linear model per cell of the feature weights of each set against the
normalised counts over all genes of the store. Cells are read and scored
in chunks, so there is no streaming switch.

## Usage

``` r
rs_lm_activity_sc(f_path, gs_list, cells_to_keep, lm_params, verbose)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- gs_list:

  List. One element per set, each a list with `idx` (integer, 0-indexed
  gene positions) and optionally `weight` (numeric, aligned with `idx`;
  absent means weight 1).

- cells_to_keep:

  Integer. Vector of indices (0-indexed) of the cells to score.

- lm_params:

  List. With `method` (`"ulm"` or `"mlm"`), `stat` (`"t"` or `"coef"`)
  and `p_values` (boolean).

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following elements

- es - Numerical matrix of cells x sets with the scores.

- p_values - Numerical matrix of cells x sets with the two-sided
  p-values, or `NULL` if not requested.

## References

Badia-i-Mompel, et al., Bioinform. Adv., 2022.
