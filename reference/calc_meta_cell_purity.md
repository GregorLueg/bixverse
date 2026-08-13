# Calculate meta cell purity

A potential metric to see how well the meta cells are aggregated is
their cell type purity. This helper function helps to plot the meta-cell
purity based on annotated cell types. These can be also just
unsupervised memberships to graph-based clustering, etc.

## Usage

``` r
calc_meta_cell_purity(
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

The `MetaCells` with added columns to the observation table:

- mc_purity - Fraction of the meta cell's cells that carry the most
  abundant label. Always added.

- mc_top_label - Name of the most abundant label. Added for
  `add_additional_info %in% c("top_label", "top_two_labels")`.

- mc_second_label - Name of the second most abundant label, `NA` for a
  pure meta cell. Added for `add_additional_info = "top_two_labels"`.

- mc_second_frac - Fraction of the meta cell's cells carrying that
  second label, `0` for a pure meta cell. Added for
  `add_additional_info = "top_two_labels"`.

- mc_entropy - Normalised Shannon entropy of the label distribution, see
  details. Added for `add_entropy = TRUE`.

## Details

Ties for the top (or second) label resolve to whichever label sorts
first, as the labels are factorised internally. The entropy is the
Shannon entropy of the label distribution within a meta cell, divided by
`log(<number of distinct labels in original_cell_type>)`, so it sits in
`[0, 1]` and stays comparable between meta cells. It is `0` if there is
only a single label in the data.
