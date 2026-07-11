# Bulk co-expression phase 2: shared helpers, unified result class, test updates

## Context

Phase 1 shipped the correctness bugs, the print methods, the NMF wiring, and a
multi-method vignette. What it deferred was the structural clean-up:

1. Every co-expression method opens with the same "grab `target_mat` or fall
   back to raw" preamble, followed by the same "guard on `detection_method`"
   block. There are five copies of each across `methods_coexp_cor.R`,
   `methods_coexp_ica.R`, `methods_coexp_dgrdl.R`, plus the new
   `methods_coexp_nmf.R`. The `&&`-vs-`||` bug in phase 1 came from exactly
   this duplication.
2. `final_results` on `BulkCoExp` has four different shapes depending on which
   method ran (data.table for Leiden, empty for CoReMo, list of matrices for
   DGRDL, list of `(S, A, ica_meta)` for ICA, list of `(modules, gene_loadings,
   sample_activity, ...)` for NMF). `get_results()` is unusable as a common
   entry point.
3. The existing tinytest files inspect `outputs$final_modules`,
   `final_results$dictionary`, `final_results$S` etc. directly. Any refactor
   of the finalisers has to bring the tests along.

Goal: extract two shared internal helpers, introduce a `BulkModuleResult` S3
class that all five method families terminate in, migrate every finaliser to
write into it, update tinytest to match, and adjust the phase-1 vignette to
use the common getters.

## Workstream 1: Shared helpers

New file `R/helpers_bulk.R` with two internals.

### `.get_bulk_target_mat(object, .verbose = TRUE)`

Replaces the five duplicated preambles:

```r
if (purrr::is_empty(S7::prop(object, "processed_data")[["processed_data"]])) {
  warning("No pre-processed data found. Defaulting to the raw data.")
  target_mat <- S7::prop(object, "raw_data")
} else {
  target_mat <- S7::prop(object, "processed_data")[["processed_data"]]
}
```

Wrap this behind one function. Suppress the warning when `.verbose = FALSE`
so the tinytest output stays clean.

### `.assert_bulk_detection_method(object, allowed, method_label)`

Replaces the five copies of:

```r
detection_method <- S7::prop(object, "params")[["detection_method"]]
if (is.null(detection_method) || !(detection_method %in% allowed)) {
  warning(sprintf(
    "This class does not seem to be set for %s module detection. Returning class as is.",
    method_label
  ))
  return(<sentinel>)
}
```

Return the resolved `detection_method` on success, `NULL` on failure. Caller
uses `if (is.null(dm)) return(object)`. Keeps the correct `||` + negation
in one place so the guard-inversion bug can't recur.

### Call sites to update

| File | Function |
|------|----------|
| `R/methods_coexp_cor.R` | `cor_module_processing`, `cor_module_tom`, `diffcor_module_processing`, `cor_module_check_epsilon`, `cor_module_graph_check_res`, `cor_module_graph_final_modules`, `cor_module_coremo_clustering`, `cor_module_coremo_stability`, `cor_module_coremo_cor_sign`, `cor_module_coremo_eigengene` |
| `R/methods_coexp_ica.R` | `ica_processing`, `ica_evaluate_comp`, `ica_stabilised_results` |
| `R/methods_coexp_dgrdl.R` | `dgrdl_grid_search`, `dgrdl_result` |
| `R/methods_coexp_nmf.R` | `nmf_bulk`, `stabilised_nmf_bulk` |

Every one of these opens with the same two blocks. The rewrite reduces each
by roughly 15 lines.

## Workstream 2: `BulkModuleResult` S3 class

New file `R/classes_bulk_module_result.R`. Constructor:

```r
new_bulk_module_result <- function(
  modules,       # data.table: gene, module_id, weight?, sign?
  factors,       # list of method-specific matrices with agnostic keys
                 #   ("gene_loadings", "sample_activity", "module_eigengenes",
                 #    "dictionary", "loadings", ...)
  method,        # character(1): one of the six detection methods
  params,        # the params_xxx() list that produced the fit
  diagnostics    # method-specific: stability DT, quality DT, laplacians,
                 # per-run losses, ...
)
```

### Methods

- `print.BulkModuleResult` — one-line dims, method, key hyperparams,
  diagnostics keys. Referenced from `print.BulkCoExp` (Section 3).
- `dim.BulkModuleResult` — `c(n_genes, n_modules)`.
- `format.BulkModuleResult` — compact one-liner for embedding in
  `print.BulkCoExp`.

### Getters (S7 generics on `S7::new_S3_class("BulkModuleResult")`)

- `get_modules(object)` — the module membership data.table.
- `get_factors(object, which = NULL)` — one factor by key, or the whole list.
- `get_diagnostics(object, which = NULL)` — same shape for diagnostics.

These are the ones referenced from `print.BulkCoExp` after the migration.

## Workstream 3: Method finaliser migration

Populate `BulkModuleResult` at the end of each method's terminal step. Keep
intermediate artefacts (epsilon sweeps, resolution grids, per-run stability)
in `outputs$*` behind their existing getters. Those are inspection surfaces,
not "the result".

| Method | Terminal step | `modules` | `factors` | `diagnostics` |
|--------|---------------|-----------|-----------|---------------|
| CoReMo | `cor_module_coremo_eigengene` | gene, module_id, sign, stability | `module_eigengenes` (sample × k), `gene_to_eigengene_cor` (gene × 1) | `stability` DT, `cluster_quality` DT |
| Leiden | `cor_module_graph_final_modules` | gene, module_id | (empty) | `resolution_used`, `modularity`, `n_subclustering_rounds` |
| Diffcor | `cor_module_graph_final_modules` (dispatches on detection method) | gene, module_id | (empty) | `resolution_used`, `n_differential_edges` |
| ICA | `ica_stabilised_results` | gene, module_id derived from `which.max(abs(S))` | `S` (gene × k gene loadings), `A` (sample × k sample activity) | `ica_meta` DT, `stability_scores` |
| DGRDL | `dgrdl_result` | gene, module_id derived from `which.max(abs(loadings))` per column | `dictionary` (sample × dict_size), `loadings` (dict_size × gene) | `feature_laplacian`, `sample_laplacian`, grid-search table |
| NMF | `nmf_bulk` / `stabilised_nmf_bulk` | already produced by `.nmf_modules_from_w` | `gene_loadings`, `sample_activity` (already present); rename existing `final_results` list keys accordingly | `final_loss`, `n_iter`, `converged`; stabilised adds `losses`, `converged`, `best_idx`, `w_all_runs`, `h_per_run` |

### ICA module derivation

Continuous loadings need a threshold policy. Reuse the phase-1 pattern from
`.nmf_modules_from_w` but on `S` (which is `k × gene` in the current code, so
we transpose): assign each gene to `which.max(abs(S_transposed))`; filter on
`min_loading` (default `0`).

### DGRDL module derivation

`loadings` is `dict_size × gene`. Assign each gene to the dictionary atom
with the largest `abs()` loading. Same `min_loading` policy.

### Print branches

`print.BulkCoExp` currently branches on `detection_method` and formats each
one directly. After the migration, the terminal step for each method
populates `final_results` with a `BulkModuleResult`. The print branches
delegate the "modules identified: N" line to `format(get_modules(x))` or
`dim(final_results)[2]` on that object. Everything upstream (preprocessing,
tested params, grid search summary) stays as it is.

## Workstream 4: Test suite updates

The existing tinytest files inspect internal shapes directly. Each one needs
a small update.

### `inst/tinytest/test_cor_modules.R`

- `expect_equal(current = final_module_data, target = expected_coremo_modules)`
  at line 288 currently reads `cor_test@outputs$final_modules`. Migrate to
  `get_modules(S7::prop(cor_test, "final_results"))`. The reference file at
  `test_data/coremo_modules.qs` will need regeneration with the new column
  shape (gene, module_id, sign, stability) — or the assertion tightened to
  the columns that actually matter.
- Graph-based Leiden block: `expected_cor_graph_res` at line 372. Assertions
  currently poke `S7::prop(cor_test, "final_results")` directly (a data.table).
  Migrate to `get_modules(S7::prop(cor_test, "final_results"))`.

### `inst/tinytest/test_ica.R`

- ICA final result assertions inspect `S7::prop(obj, "final_results")$S` and
  `$A` directly. Migrate to `get_factors(..., which = "S")` /
  `get_factors(..., which = "A")`. Add a module-membership assertion using
  the new `get_modules()`.

### `inst/tinytest/test_sparse_dict_learn.R`

- DGRDL assertions inspect `final_results$dictionary` and `$loadings`. Same
  migration to `get_factors()`. Add a `get_modules()` assertion.

### New test: `inst/tinytest/test_bulk_nmf.R`

Currently absent. Cover:
- `nmf_bulk` on `synthetic_bulk_cor_matrix()` output, `k = 5`. Assert
  `is(get_modules(...), "data.table")`, expected loadings shape, non-negative
  entries in `gene_loadings` and `sample_activity`.
- `stabilised_nmf_bulk` with `n_runs = 3`. Assert diagnostics has `losses`,
  `converged`, `best_idx`.
- Recall vs planted modules ≥ 0.5 (soft floor — synthetic data is easy).

### New test: `inst/tinytest/test_ica_cv_branch.R`

The `cross_validate = TRUE` branch of `ica_stabilised_results` was crashing
on `X_raw` until phase 1. Lock that fix in with a minimal test that runs the
CV branch on a small synthetic and asserts it returns a populated
`BulkModuleResult` without error.

### Reference data files under `inst/tinytest/test_data/`

- `coremo_modules.qs` — regenerate as the modules data.table (gene,
  module_id, sign, stability) rather than the current wide shape.
- `cor_graph_final_res.qs` — regenerate as the modules data.table.

Regeneration is a one-liner in each test file: rerun the pipeline, then
`qs2::qs_save(get_modules(final), "./test_data/coremo_modules.qs")`. Do that
inside the test file behind a `if (identical(Sys.getenv("REGEN"), "1"))`
guard so it doesn't run on every invocation.

## Workstream 5: Vignette adjustment

`vignettes/bulk_coexpression_modules.qmd` was written against the phase-1
per-method final shapes. After the migration:

- CoReMo section: replace `get_outputs(coexp_coremo)$final_modules` with
  `get_modules(S7::prop(coexp_coremo, "final_results"))`.
- Leiden section: `get_results(coexp_leiden)` still works (inherited generic)
  but now returns a `BulkModuleResult`. Update to `get_modules(...)`.
- ICA section: replace `S7::prop(obj, "final_results")$S` with
  `get_factors(..., which = "S")` and derive modules via `get_modules(...)`.
- DGRDL section: same pattern for `dictionary` / `loadings`.
- NMF section: the phase-1 `get_nmf_modules()` / `get_nmf_gene_loadings()`
  wrappers stay as thin shims that call `get_modules()` /
  `get_factors(..., "gene_loadings")` under the hood. Keeps the public API
  stable.

The recovery-plot section at the end reduces to one loop over methods, all
using the same `get_modules()` call. That is the visible payoff.

## Order of implementation

1. **Helpers first** (Workstream 1). Small, independent, unlocks the guard
   uniformity.
2. **`BulkModuleResult` class + methods** (Workstream 2). Standalone; nothing
   in existing code touches it yet.
3. **Migrate one finaliser end-to-end** (start with CoReMo — smallest test
   footprint). Update the two `test_data/*.qs` files. Confirm tests pass.
4. **Migrate remaining finalisers** in this order: Leiden, ICA, DGRDL, NMF.
   After each, run the corresponding tinytest file.
5. **Add new tests** (`test_bulk_nmf.R`, `test_ica_cv_branch.R`).
6. **Vignette pass** last, so the getters it calls all exist and are stable.
7. **Full `tinytest::run_test_dir()` + `devtools::check()`** to close out.

## Files touched

**New:**
- `R/helpers_bulk.R`
- `R/classes_bulk_module_result.R`
- `inst/tinytest/test_bulk_nmf.R`
- `inst/tinytest/test_ica_cv_branch.R`

**Modified:**
- `R/methods_coexp_cor.R` — swap preambles for helpers; rewrite the two
  terminal steps (`cor_module_graph_final_modules`, `cor_module_coremo_eigengene`) to build `BulkModuleResult`.
- `R/methods_coexp_ica.R` — swap preambles; rewrite `ica_stabilised_results`
  to build `BulkModuleResult`.
- `R/methods_coexp_dgrdl.R` — swap preambles; rewrite `dgrdl_result` to
  build `BulkModuleResult`.
- `R/methods_coexp_nmf.R` — swap preambles; rewrite `nmf_bulk` and
  `stabilised_nmf_bulk` finalisers; keep the `get_nmf_*` wrappers as shims
  over `get_modules` / `get_factors`.
- `R/classes_bulk.R` — `print.BulkCoExp` branches delegate module counts to
  `format(get_modules(x))`; `get_results()` docs note the new return type.
- `inst/tinytest/test_cor_modules.R`, `test_ica.R`, `test_sparse_dict_learn.R`
  — assertions migrated to the shared getters.
- `inst/tinytest/test_data/coremo_modules.qs`,
  `test_data/cor_graph_final_res.qs` — regenerated behind `REGEN=1` guards.
- `vignettes/bulk_coexpression_modules.qmd` — assertions migrated to
  `get_modules()` / `get_factors()`.

## Verification

1. `rextendr::document()` clean (no rust changes expected in phase 2).
2. `devtools::document()` clean.
3. `devtools::check()` — no NOTEs about undefined globals, no warnings from
   the migrated preambles.
4. `tinytest::run_test_dir("inst/tinytest")` — every file passes. The
   pre-existing `test_cistarget.R` sourcing issue that surfaced during
   phase 1 is a separate concern (env not seeing `run_cistarget`) and is out
   of scope here.
5. Render the phase-1 vignette with `quarto render vignettes/bulk_coexpression_modules.qmd`;
   confirm every chunk still runs and the recovery bar plot still shows five
   method scores.
6. Interactive: build a `BulkCoExp`, run each method, call `get_modules`,
   `get_factors`, `get_diagnostics`. Confirm consistent return shapes.

## Migration notes

- `BulkModuleResult` is S3 by choice. R6 would give inheritance we don't
  need; S7 would tie the result to the S7 class registry when the natural
  home for a result-of-a-fit is a plain named list with a class. S3 fits.
- The `get_nmf_*` getters from phase 1 stay as thin wrappers so existing
  user code doesn't break. Same policy applies to any user code we don't
  own that pokes at `final_results$dictionary` and similar — those users
  can migrate to `get_factors()` on their own timeline.
- No CRAN-visible API break. The visible change is that `get_results()`
  returns a `BulkModuleResult` instead of raw shapes; that's an addition,
  not a removal.

## Not in scope for phase 2

- Retrofitting `BulkModuleResult` onto contrastive PCA. It doesn't produce
  modules in the same sense (the output is factors + loadings, not gene
  groupings). Leave `contrastive_pca` alone.
- Deprecating the per-method getters (`get_epsilon_res`,
  `get_resolution_res`, `get_grid_search_res`, `get_ica_stability_res`).
  Those cover intermediate artefacts, not the terminal result, and they
  stay.
- Refactoring `BulkDge`. Its `final_results` is a different beast and is
  fine as is.
