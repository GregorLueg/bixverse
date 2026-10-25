# Generate meta cells via MetaCells2 (divide and conquer)

**\[experimental\]** Rust port of the metacells divide and conquer
pipeline. Works on the raw counts directly, so no embedding or kNN graph
is needed. Cells are split into piles, each pile gets its own feature
selection, similarity graph and partition, and outliers are re-run in
further rounds. Cells sharing a rare gene module are piled together
first unless switched off in `dac_params`.

## Usage

``` r
rs_metacells2(
  f_path,
  cells_to_keep,
  mc2_params,
  dac_params,
  target_size,
  seed,
  verbose
)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- cells_to_keep:

  Integer vector. 0-indexed positions in the count file of the cells to
  include.

- mc2_params:

  List. The per-pile parameters, see
  [`params_sc_metacells2()`](https://gregorlueg.github.io/bixverse/reference/params_sc_metacells2.md).

- dac_params:

  List. The divide and conquer parameters, see
  [`params_sc_metacells2_dac()`](https://gregorlueg.github.io/bixverse/reference/params_sc_metacells2_dac.md).

- target_size:

  Numeric. Library size to normalise the aggregated counts to.

- seed:

  Integer. For reproducibility.

- verbose:

  Boolean. Controls verbosity.

## Value

A list with the following elements

- assignments - The meta cell assignments in the original cell space.

- aggregated - The aggregated raw and normalised counts in CSR form.

- dissolved - Logical per included cell. Was the cell's candidate
  dissolved in its last pile.

- level - Integer per included cell. Outlier round that placed it.

- rare_module_of_cell - Integer per included cell. 0-indexed rare gene
  module, `-1` for none.

- rare_module_of_gene - Integer per gene. 0-indexed rare gene module,
  `-1` for none.

- selected_genes - Integer. 1-indexed genes selected in any pile plus
  the rare module genes.

## References

Ben-Kiki, et al., Genome Biol., 2022.
