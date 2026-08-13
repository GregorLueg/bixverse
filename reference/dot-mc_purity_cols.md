# Build the meta cell purity columns

Shared body of
[`calc_meta_cell_purity()`](https://gregorlueg.github.io/bixverse/reference/calc_meta_cell_purity.md)
and
[`get_meta_cell_purity()`](https://gregorlueg.github.io/bixverse/reference/get_meta_cell_purity.md).
Resolves the label vector against the right index space and tabulates
the label distribution per meta cell.

## Usage

``` r
.mc_purity_cols(object, original_cell_type, add_additional_info, add_entropy)
```

## Arguments

- object:

  `MetaCells` class.

- original_cell_type:

  Character vector. The original cell type annotations, either in full
  obs order or in QC-passing order. Which one was passed is inferred
  from the length.

- add_additional_info:

  String. One of `c("none", "top_label", "top_two_labels")`. Already
  resolved by [`match.arg()`](https://rdrr.io/r/base/match.arg.html) in
  the caller.

- add_entropy:

  Boolean. Add the normalised Shannon entropy.

## Value

A named list of atomic vectors, one element per meta cell in `obs_table`
row order:

- mc_purity - Fraction of the meta cell's cells that carry the most
  abundant label. Always present.

- mc_top_label - Name of the most abundant label. Present for
  `add_additional_info %in% c("top_label", "top_two_labels")`.

- mc_second_label - Name of the second most abundant label, `NA` for a
  pure meta cell. Present for `add_additional_info = "top_two_labels"`.

- mc_second_frac - Fraction of the meta cell's cells carrying that
  second label, `0` for a pure meta cell. Present for
  `add_additional_info = "top_two_labels"`.

- mc_entropy - Normalised Shannon entropy of the label distribution.
  Present for `add_entropy = TRUE`.
