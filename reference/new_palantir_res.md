# Helper function to generate the Palantir results

Takes the raw Rust output of
[`rs_palantir()`](https://gregorlueg.github.io/bixverse/reference/rs_palantir.md)
and maps every cell index back onto the cell names the kNN graph was
built over.

## Usage

``` r
new_palantir_res(rs_res, used_cells, modality, terminal_labels = NULL)
```

## Arguments

- rs_res:

  List. The raw return of
  [`rs_palantir()`](https://gregorlueg.github.io/bixverse/reference/rs_palantir.md).

- used_cells:

  Character vector. The cells the kNN graph was generated over, in kNN
  row order.

- modality:

  String. The modality the kNN graph came from.

- terminal_labels:

  Optional named character vector. The terminal states as the caller
  supplied them, i.e. lineage label to cell name. When given, the labels
  become the column names of `branch_probs` and the names of
  `terminal_states`. Detected terminal states have no labels to carry,
  so this is `NULL` for those runs.

## Value

Generates the `PalantirRes` class.
