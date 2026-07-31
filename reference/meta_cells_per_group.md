# Generate source-pure meta cells and merge them

Splits `object` by `group_col`, optionally runs `pipeline` on each
subset, generates meta cells per group and merges the results into a
single
[`MetaCells()`](https://gregorlueg.github.io/bixverse/reference/MetaCells.md)
object. This gives you meta cells that never mix cells from two
patients/samples, while still being one object you can run SCENIC,
AUCell or NMF over.

The meta cell generators need an embedding, so `pipeline` will normally
be `step_hvg_sc() %>>% step_pca_sc() %>>% step_neighbours_sc()` unless
every subset already carries one.

## Usage

``` r
meta_cells_per_group(
  object,
  group_col,
  method = c("bootstrapped", "seacells", "supercells"),
  mc_params = list(),
  pipeline = NULL,
  groups = NULL,
  feature_space = c("intersect", "union"),
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`.

- group_col:

  String. Column in obs used to split.

- method:

  String. One of `c("bootstrapped", "seacells", "supercells")`. Picks
  the meta cell generator.

- mc_params:

  Named list. Arguments passed on to the generator, e.g.
  `list(sc_meta_cell_params = params_sc_bt_metacells(), target_size = 1e5)`.

- pipeline:

  Optional `ScPipeline` applied to each subset before the meta cells are
  generated.

- groups:

  Optional character vector. Restrict to these group values; if `NULL`,
  all unique values of `group_col` are used.

- feature_space:

  String. One of `c("intersect", "union")`. Passed to
  [`merge_meta_cells()`](https://gregorlueg.github.io/bixverse/reference/merge_meta_cells.md).
  Irrelevant here as all groups share the gene space of the parent
  object.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A merged
[`MetaCells()`](https://gregorlueg.github.io/bixverse/reference/MetaCells.md)
object with a `source_id` column in its observation table.

## Examples

``` r
if (FALSE) { # \dontrun{
prep <- step_hvg_sc() %>>% step_pca_sc(no_pcs = 20L) %>>% step_neighbours_sc()
mc <- meta_cells_per_group(
  object = sc_obj,
  group_col = "patient_id",
  method = "bootstrapped",
  pipeline = prep
)
} # }
```
