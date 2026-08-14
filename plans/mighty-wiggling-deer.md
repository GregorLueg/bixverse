# Fix the deserialised-Matrix parsing bug, test the meta cell merge, de-duplicate the sc test fixtures

## Context

Pushing a `MetaCells` object through Rust after a `qs2` round trip fails with:

```
R list parsing failed: nrow missing or not a non-negative whole number
```

The suspicion was `merge_meta_cells()` and an integer/double mismatch on `nrow`.
Neither holds up.

`Matrix` is listed in DESCRIPTION `Imports:` but nothing in `NAMESPACE` imports
from it, so `library(bixverse)` never loads the Matrix namespace. When a
`dgRMatrix` is deserialised into such a session, the **first** S4-dispatched
call on it silently misfires while R lazily loads the namespace. Reproduced on
the installed 0.4.8 build:

```r
# fresh session
library(bixverse)
m <- qs2::qs_read("mc.qs2")
is.null(nrow(m@data$raw))                 # TRUE  <- the bug
m@data$raw[1:5, , drop = FALSE]           # object of type 'S4' is not subsettable
find_hvg_sc(m, hvg_no = 5L)               # R list parsing failed: nrow missing ...
merge_meta_cells(list(a, b))              # Result must be length 1, not 0
```

`mc_counts_to_list()` then hands `nrow = NULL` to Rust, and
`list_to_sparse_matrix()` reports it as missing. A second call succeeds, because
the first one triggered the namespace load. Not merge-specific: an unmerged
`MetaCells` fails the same way. `saveRDS`/`readRDS` behaves identically, so this
is not a qs2 quirk.

`as.numeric()` around `nrow`/`ncol` fixes nothing. `bixverse-rs` 0.4.3 (the
locked version, `r_rust_interface.rs:653`) already accepts an integer or a whole
double, and `as.numeric(NULL)` is `numeric(0)`, still not a scalar.

Separately, `merge_meta_cells()` and its helpers have zero test coverage, and the
single cell test files carry five blocks of copy-pasted fixture scaffolding
across ~22 files.

## Part 1: the fix

**`R/bixverse-package.R`** - add a real import so the Matrix namespace is loaded
when bixverse is:

```r
## usethis namespace: start
#' @importFrom lifecycle deprecated
#' @importFrom Matrix sparseMatrix
## usethis namespace: end
NULL
```

`sparseMatrix` is the natural pick: `get_meta_cell_matrices()`
(`R/functions_single_cell.R:218`) already calls it. Then `devtools::document()`
so `NAMESPACE` picks up `importFrom(Matrix,sparseMatrix)`.

This is the whole fix. It covers every bare `nrow()` / `ncol()` / `dim()` / `[` /
`t()` on a Matrix object in the package, not just the one in
`mc_counts_to_list()`.

**`R/classes_meta_cell.R:428-429`** - while here, read the dimensions off the
slot rather than through dispatch, matching what `sparse_mat_to_list()`
(`R/helpers.R:575`) already does:

```r
    nrow = x@Dim[1],
    ncol = x@Dim[2]
```

Hardening, not the fix. It would not have saved the subsetting path or
`.merge_mc_counts()`.

Do **not** add `as.numeric()`.

## Part 2: shared test fixtures

Fill the empty `inst/tinytest/helper_sc.R` (committed in `80a1c785`, 0 bytes,
never sourced). tinytest has no helper convention, but `run_test_file()`
`setwd()`s to the test directory and evaluates each top-level expression in one
env, so the first line of each test file works:

```r
source("helper_sc.R", local = TRUE)
```

Verify this under both `tinytest::run_test_file()` and
`tinytest::test_package("bixverse")` before migrating anything - the latter runs
against the installed `<pkg>/tinytest/` directory.

Helpers to write, all with roxygen-free but commented signatures and
`checkmate` validation on their arguments (they are test code, keep them short):

| Helper | Replaces | Notes |
|---|---|---|
| `sc_test_dir(name)` | Block A (22 files) + `get_obj_dir()` (3 verbatim copies at `test_sc_obj_manipulation.R:22`, `test_sc_references.R:22`, `test_sc_sctype.R:21`) | `file.path(tempdir(), name)`, `dir.create(recursive = TRUE)`, `stopifnot(dir.exists(...))`, returns the path. One function covers both the top-level dir and the nested per-object dirs. |
| `sc_test_fixture()` | Block B (10 files) | Returns a list: thresholds (`min_lib_size`, `min_genes_exp`, `min_cells_exp`, `hvg_to_keep`, `no_pcs`), `data` from `generate_single_cell_test_data()`, plus `genes_pass` / `cells_pass`. **No `expect_*` inside** - tinytest only registers expectations returned at top level. Keep the two sanity assertions in one test file (`test_sc_processing.R`), drop the five other copies. |
| `sc_test_object(dir, fixture)` | Block C (27 call sites, 19 files) | `SingleCells(dir_data = dir)` + `params_sc_min_quality(...)` + `load_r_data(..., streaming = 0L, .verbose = FALSE)`. |
| `sc_test_prepped(object, fixture)` | the `find_hvg_sc` -> `calculate_pca_sc` -> `find_neighbours_sc` chain (~8 files) | Optional convenience on top of `sc_test_object()`. |
| `sc_batch_fixture(strength)` | Block E (`test_sc_batch_corr.R:25-30`, `test_sc_subset_batch_corr.R:17-22`) | The 4-cell-type `cell_markers` list + `generate_single_cell_test_data(params_sc_synthetic_data(n_cells = 900L, n_batches = 3L, ...))`. |
| `sc_test_cleanup(...)` | the trailing `on.exit(unlink(...))` idiom | Plain `unlink(..., recursive = TRUE, force = TRUE)`. The current top-level `on.exit()` does not actually defer (tinytest evals expressions one at a time, so it fires immediately) and never registers if a test above it errors. A plain call as the last statement is honest about what it does. |

Then migrate the sc/mc test files to use them. Mechanical; the generator itself
(`generate_single_cell_test_data()`, `R/data_single_cell.R:35`) already lives in
`R/` and stays there.

Files that deviate and should be left alone: `test_sc_lr.R` (pure data.table
fixtures), the first ~400 lines of `test_sc_trajectory.R` (hand-rolled
manifolds), `test_sc_doublet_detection.R` (builds doublets by summing count
rows).

Five files currently leak their temp dir for the session
(`test_sc_cache_state.R`, `test_sc_obj_manipulation.R`, `test_sc_plot_helpers.R`,
`test_sc_sctype.R`, `test_sc_trajectory.R`) - give them the cleanup call too.

## Part 3: tests for the meta cell merge

New file `inst/tinytest/test_mc_merge.R`. Two `MetaCells` objects built through
`sc_test_object()` -> `SingleCellsSubset` -> `generate_bt_meta_cells_sc()`, plus a
small hand-built pair for the gene space edge cases (constructor takes
`meta_cell_data` / `var_data` / `meta_cell_method` directly, see
`R/classes_meta_cell.R:65`).

Cover:

- **Gene space.** `feature_space = "intersect"` and `"union"` over inputs with
  differing gene sets *and* differing gene order, so the non-identity remap and
  the `is.unsorted()` re-sort branch in `.merge_mc_counts()`
  (`R/functions_meta_cell.R:564-589`) both fire. Assert the merged matrix equals
  a dense `rbind` reference built with `match(target_genes, gene_ids)`, and that
  `validObject()` passes on both `data$raw` and `data$norm`. I verified both
  branches are numerically correct today, so this is a regression lock.
- **Structural zeros.** Under `"union"`, genes absent from an input stay zero for
  its meta cells.
- **obs table.** `source_id` present, `meta_cell_id` prefixed as
  `<source>__<id>`, `meta_cell_idx` renumbered `1:.N`, column order
  `meta_cell_idx, meta_cell_id, source_id`, `original_cell_idx` carried over
  untouched.
- **Provenance.** `is_merged` TRUE, `other_data$sources`,
  `original_assignment$per_source` holding only the three scalars,
  `assignments = NULL`, `cells_to_keep = NULL`, `dims` matching the merged
  matrices. Shared-parent branch (`R/functions_meta_cell.R:127-141`): equal
  `n_cells` across sources gives one `n_cells` and a bitmap-derived
  `n_unassigned`, unequal sums both.
- **Caches dropped.** `sc_cache` reset - PCA/kNN set on the inputs is gone.
- **Errors and warnings.** Mismatched `meta_cell_method`; empty intersection;
  non-unique `source_ids`; `prefix_ids = FALSE` with duplicated ids; the
  dropped-genes warning; the non-shared obs column warning; the `other_data`
  warning.
- **Downstream.** `find_hvg_sc()`, `calculate_pca_sc()`, `find_neighbours_sc()`,
  `aucell_sc()`, `scenic_grn_sc()`, `nmf_sc()` all run on the merged object and
  return the expected shapes (`dims[1]` rows). This is the "does it survive the
  Rust boundary" check.
- **Refusals.** `calc_diffusion_coordinates()` and `calc_manifold_metrics()`
  error on a merged object (`.mc_artefact_rows()`, `R/functions_meta_cell.R:208`).
- **Pipeline.** One `meta_cells_per_group()` end-to-end case in
  `test_sc_pipelines.R` or the new file, since it is the only in-package caller
  (`R/pipelines_sc.R:412`).

Check each generic's signature in `R/base_generics_sc.R` before writing the
call - `sc_pipeline()` takes no arguments and steps are chained with `%>>%`, not
passed positionally.

## Part 4: regression test for the namespace bug

The lazy-load failure cannot be reproduced inside a running test process, since
Matrix is loaded by then. Lock the fix directly instead, in `test_utils.R` or the
new file:

```r
expect_true(
  current = "Matrix" %in% names(getNamespaceImports("bixverse")),
  info = "bixverse imports from Matrix so the namespace loads eagerly"
)
```

Optionally, one `at_home()`-free subprocess check that a saved `MetaCells`
survives a fresh session, guarded on `requireNamespace("qs2")` (already in
Suggests) and on `Sys.which("Rscript")`. Skip it if it turns out flaky on CI -
the import assertion is the load-bearing one.

## Verification

```sh
DEV_BUILD=true R CMD INSTALL --preclean .
air format .
```

```r
devtools::document()                                  # NAMESPACE gains Matrix
tinytest::run_test_file("inst/tinytest/test_mc_merge.R")
tinytest::test_package("bixverse")                    # full suite, migrated helpers
```

End-to-end check of the original failure, in a **fresh** process:

```sh
Rscript -e 'library(bixverse); m <- qs2::qs_read("mc.qs2"); find_hvg_sc(m, hvg_no = 5L)'
Rscript -e 'library(bixverse); a <- qs2::qs_read("mc_single.qs2"); merge_meta_cells(list(x = a, y = a))'
```

Both currently fail on the first call and pass on a retry. After the import they
must pass first time.
