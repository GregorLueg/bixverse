# Calculate meta cell purity without mutating object state

Like
[`calc_meta_cell_purity()`](https://gregorlueg.github.io/bixverse/reference/calc_meta_cell_purity.md)
but does not mutate `object`. Returns a data.table with the purity per
meta cell instead of stamping the columns onto the observation table.
Useful for comparing several label columns, or for sweeping meta cell
resolutions, without permanently annotating the object.

## Usage

``` r
get_meta_cell_purity(
  object,
  original_cell_type,
  add_additional_info = c("none", "top_label", "top_two_labels"),
  add_entropy = FALSE
)
```

## Arguments

- object:

  `MetaCells` class.

- original_cell_type:

  Character vector. The original cell type annotations of the object the
  meta cells came from. Either in the row order of its full (unfiltered)
  obs table, i.e. `get_sc_obs(x)$<column>`, or of the QC-passing cells
  only, i.e. `get_sc_obs(x, filtered = TRUE)$<column>`. Which one you
  passed is inferred from the length, so a vector matching neither is an
  error rather than a silently wrong purity.

- add_additional_info:

  String. Which label information to add on top of the purity. One of
  `c("none", "top_label", "top_two_labels")`. Defaults to `"none"`, i.e.
  the purity only.

- add_entropy:

  Boolean. Shall the normalised Shannon entropy of the label
  distribution be added as a diversity measure. Defaults to `FALSE`.

## Value

A data.table in observation table row order:

- meta_cell_idx - Index of the meta cell. Always returned.

- meta_cell_id - Identifier of the meta cell. Always returned.

- mc_purity - Fraction of the meta cell's cells that carry the most
  abundant label. Always returned.

- mc_top_label - Name of the most abundant label. Returned for
  `add_additional_info %in% c("top_label", "top_two_labels")`.

- mc_second_label - Name of the second most abundant label, `NA` for a
  pure meta cell. Returned for `add_additional_info = "top_two_labels"`.

- mc_second_frac - Fraction of the meta cell's cells carrying that
  second label, `0` for a pure meta cell. Returned for
  `add_additional_info = "top_two_labels"`.

- mc_entropy - Normalised Shannon entropy of the label distribution, see
  details. Returned for `add_entropy = TRUE`.

## Details

Ties for the top (or second) label resolve to whichever label sorts
first, as the labels are factorised internally. The entropy is the
Shannon entropy of the label distribution within a meta cell, divided by
`log(<number of distinct labels in original_cell_type>)`, so it sits in
`[0, 1]` and stays comparable between meta cells. It is `0` if there is
only a single label in the data.
