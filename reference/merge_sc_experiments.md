# Merge multiple `SingleCells` experiments into one

Merges N existing `SingleCells` objects into a freshly constructed
target object. The feature space of the result is the **intersection**
of the input gene sets. Each input's `cells_to_keep` filter is honoured
(i.e. only cells with `to_keep = TRUE` in the input's obs are carried
over).

If `renormalise = FALSE`, the stored `data_norm` values are copied
through unchanged. This is valid only when all inputs were normalised
against the same `target_size`. If the gene intersection is much smaller
than the individual input gene sets, the inherited `data_norm` becomes a
lossy approximation (it was computed against the pre-intersection
library size). In that case set `renormalise = TRUE` to recompute
`data_norm` against the surviving raw counts using
`sc_qc_param$target_size`.

Obs columns are intersected across inputs. The result obs gains an
`exp_id` column. Inputs that already have an `exp_id` column are
rejected. The `sc_cache` and `sc_map` of the target are populated fresh;
any PCA, kNN, sNN or HVG state on the inputs is not carried over and
must be re-run.

## Usage

``` r
merge_sc_experiments(
  target,
  inputs,
  exp_ids,
  renormalise = FALSE,
  sc_qc_param = params_sc_min_quality(),
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
)
```

## Arguments

- target:

  A freshly constructed `SingleCells` pointing at the output directory.

- inputs:

  List of `SingleCells` objects to merge. Length \>= 2.

- exp_ids:

  Character vector of experiment identifiers, one per input. Must be
  unique.

- renormalise:

  Boolean. Whether to recompute `data_norm` against
  `sc_qc_param$target_size`. Defaults to `FALSE`.

- sc_qc_param:

  List. Output of
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).
  Only `target_size` is consulted here; no QC filtering is applied
  during merge.

- csc_mem_gb:

  Optional numeric. Memory in GB for the buffers of the cell-to-gene
  (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL` (default)
  converts in one pass and holds the whole matrix. Set a cap for large
  data sets; every extra phase re-reads the cell file once.

- streaming, batch_size, max_genes_in_memory, cell_batch_size:

  Replaced by `csc_mem_gb` and ignored. **\[deprecated\]**

- .verbose:

  Boolean.

## Value

The populated target `SingleCells`.

## Examples

``` r
# \donttest{
# two synthetic experiments merged over their shared gene space
sc_a <- demo_single_cells(prepped = FALSE, seed = 1L)
sc_b <- demo_single_cells(prepped = FALSE, seed = 2L)
merged_dir <- tempfile("bixverse_merged")
dir.create(merged_dir)

merged <- merge_sc_experiments(
  target = SingleCells(dir_data = merged_dir),
  inputs = list(sc_a, sc_b),
  exp_ids = c("exp_a", "exp_b"),
  .verbose = FALSE
)
table(unlist(merged[["exp_id"]]))
#> 
#> exp_a exp_b 
#>   500   500 

unlink(
  c(sc_a@dir_data, sc_b@dir_data, merged_dir),
  recursive = TRUE,
  force = TRUE
)
# }
```
