# Helper function to generate the Bonsai results

Takes the raw Rust output of
[`rs_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_bonsai.md)
and maps the leaves back onto the cell names and the genes back onto
their identifiers.

## Usage

``` r
new_bonsai_tree(
  rs_res,
  cell_idx,
  cell_names,
  genes_in,
  gene_ids,
  bonsai_params,
  leaf_sizes = NULL
)
```

## Arguments

- rs_res:

  List. The raw return of
  [`rs_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_bonsai.md)
  or
  [`rs_mc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/rs_mc_bonsai.md).

- cell_idx:

  Integer. The cells or metacells the tree was built over (0-indexed!),
  in leaf order.

- cell_names:

  Character vector. The names of those cells or metacells.

- genes_in:

  Integer. The genes that went in (0-indexed!).

- gene_ids:

  Character vector. Identifiers of `genes_in`, same order.

- bonsai_params:

  List. The parameters of the run, see
  [`params_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/params_sc_bonsai.md).

- leaf_sizes:

  Optional integer. Cells behind each leaf. `NULL` means one each, i.e.
  single cells.

## Value

Generates the `BonsaiTree` class.
