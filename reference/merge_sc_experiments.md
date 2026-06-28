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
  streaming = 1L,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
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

- streaming:

  Integer. `0` -\> no streaming, `1` -\> light streaming, `2` -\> heavy
  streaming with memory upper boundaries. This enables you to control
  the memory pressure during ingestion.

- batch_size:

  Integer. Batch size when `streaming = 1L`.

- max_genes_in_memory:

  Integer. How many genes shall be held in memory at a given point.
  Defaults to `2000L`. Only relevant if streaming is set to `2`.

- cell_batch_size:

  Integer. How big are the batch sizes for the cells in the
  transformation from the cell-based to gene-based format. Defaults to
  `100000L`. Only relevant if streaming is set to `2`.

- .verbose:

  Boolean.

## Value

The populated target `SingleCells`.
