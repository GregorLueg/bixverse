# Wire one-vs-many AUROC markers into bixverse

## Context

`bixverse-rs` gained `calculate_dge_one_vs_many_auroc` on the local, unpushed
branch `test-coverage` (commit `8c73ad8`). It scores one reference cell group
against each rival group separately and summarises per gene across rivals. That
is the specific-marker question, and the R package cannot ask it today.

What exists on the R side is `find_markers_sc` (two groups) and
`find_all_markers_sc` (group vs. everything else pooled). The pooled version is
dominated by whichever rival contributes the most cells, so a gene that is high
in the reference and equally high in one small rival still looks like a clean
marker. The new function fixes that: `min_auroc`, `median_auroc` and `min_rank`
across rivals, plus Simes and max-p combinations, i.e. scran's `scoreMarkers`
summaries.

Upstream is purely additive (one public fn, one result struct, four error
variants), so nothing existing changes behaviour.

## Upstream signature

```rust
pub fn calculate_dge_one_vs_many_auroc<S: SingleCellReading>(
    reader: &S,
    ref_indices: &[usize],       // 0-based cell indices
    other_indices: &[Vec<usize>],// 0-based, one vector per rival
    min_proportion: f32,
    alternative: &str,           // "twosided" | "greater" | "less"
    verbose: usize,
) -> Result<DgeAurocMultiRes, BixverseErrors>
```

`DgeAurocMultiRes` (`src/single_cell/sc_analysis/dge_pathway_scores.rs:63`):

- Per rival, each `Vec<Vec<_>>` of length `n_rivals`, inner length `n_kept`:
  `auroc`, `lfc`, `prop_other`, `z_scores`, `p_vals`, `fdr`
- Per gene, length `n_kept`: `prop_ref`, `median_auroc`, `min_auroc`,
  `mean_auroc`, `max_auroc`, `worst_comparison` (0-based rival index),
  `min_rank`, `simes_p`, `simes_fdr`, `max_p`, `max_p_fdr`
- `genes_to_keep`: logical over **all** genes

Gene filtering is global across all groups, so every arm's FDR is computed over
the same gene set. Degenerate case: no gene clears `min_proportion`, then
`genes_to_keep` is all `FALSE` and every result vector is empty.

## Step 1: point the Cargo patch at the local checkout

`src/rust/Cargo.toml:42` currently reads:

```toml
[patch.crates-io]
bixverse-rs = { git = "https://github.com/GregorLueg/bixverse-rs", branch = "release-0.4.3" }
```

`test-coverage` is local only (no `origin/test-coverage`), so the git patch
cannot see the new function. Swap to a path override. This is a git worktree, so
the relative path is six levels up, not the usual two:

```toml
[patch.crates-io]
bixverse-rs = { path = "../../../../../../bixverse-rs" }
```

This edit is local scaffolding. Revert it to the git/version form before
committing, and bump the crate version + `[dependencies]` version together once
`bixverse-rs` is released.

## Step 2: extendr binding

`src/rust/src/single_cell/r_sc_analysis.rs`, in the existing `// DGEs` section
directly under `rs_calculate_dge_mann_whitney`. No new imports:
`dge_pathway_scores::*` is already glob-imported at line 6.

```rust
fn rs_calculate_dge_one_vs_many(
    f_path: String,
    cell_indices_ref: &[i32],
    cell_indices_other: List,
    min_prop: f64,
    alternative: String,
    verbose: usize,
) -> extendr_api::Result<List>
```

Register as `fn rs_calculate_dge_one_vs_many;` in the `extendr_module!` block
under the `// dge` comment (line 33).

Body follows `rs_calculate_dge_mann_whitney` verbatim for the reader
construction (`ParallelSparseReader::new(&f_path).to_extendr()?`) and error
handling (`.to_extendr()?`). For `cell_indices_other`, reuse the
`List` -> `Vec<Vec<usize>>` loop already used for `gs_list` in
`rs_module_scoring` (`r_sc_analysis.rs:207`).

Return a flat named list in two blocks. extendr wants `f64`/`i32`, so `f32` and
`usize` need converting the same way the pairwise binding does.

Per-rival, length `n_rivals * n_kept`, rival-major (flatten the outer `Vec`
first so the ordering is unambiguous):

- `comparison` (`i32`, **0-based** rival index, repeated `n_kept` times each)
- `auroc`, `lfc`, `prop_other` (`f32` -> `f64`)
- `z_scores`, `p_values`, `fdr` (already `f64`)

Per gene, length `n_kept`: `prop_ref`, `median_auroc`, `min_auroc`,
`mean_auroc`, `max_auroc` (`f32` -> `f64`), `worst_comparison` and `min_rank`
(`usize` -> `i32`), `simes_p`, `simes_fdr`, `max_p`, `max_p_fdr`.

Plus `genes_to_keep` (logical, length = total genes).

Index convention: `comparison` and `worst_comparison` stay **0-based**, matching
how `rs_make_milor_nhoods` returns `nhoods_i`/`nhoods_j`. The R side adds 1.
Document this in the roxygen `@param`/`@return` block, which follows the
existing style: `` `r lifecycle::badge("experimental")` ``, `@export`,
`@keywords internal`.

## Step 3: no new parameter wrapper

Three scalars (`method`, `alternative`, `min_prop`) do not earn a `params_*()`
constructor, and both DGE siblings pass them flat. Nothing to add to
`R/param_wrappers.R` or `R/checkmate_extensions.R`.

## Step 4: S3 result class

`R/classes_single_cell_others.R`, new subsection. Copy the shape of
`new_sc_type_cell_results` / `print.ScTypeCellResults` (line 2614 onwards).

```r
new_sc_specific_markers(summary, per_comparison, params)  # @keywords internal
print.ScSpecificMarkers
```

`summary`: `ref_grp`, `gene_id`, `prop_ref`, `median_auroc`, `min_auroc`,
`mean_auroc`, `max_auroc`, `worst_rival`, `min_rank`, `simes_p`, `simes_fdr`,
`max_p`, `max_p_fdr`.

`per_comparison`: `ref_grp`, `rival_grp`, `gene_id`, `auroc`, `lfc`, `prop_ref`,
`prop_rival`, `z_scores`, `p_values`, `fdr`.

`params`: the resolved `method`, `alternative`, `min_prop`,
`column_of_interest`, `downsampling`, `seed`.

`worst_rival` and `rival_grp` are resolved to group names in R (`+ 1` on the
Rust index, then index into the rival name vector for that arm). No getters
beyond `$summary` / `$per_comparison`.

Register `"ScSpecificMarkers"` in the print-method vector in `R/zzz.R`.

## Step 5: S7 generic and method

`R/methods_sc_dge_diffabund.R`, new `### find specific markers` subsection after
the `find_all_markers_sc` block (ends line 293). Dispatch on `ScOrScSubset`
(`R/classes_union.R:10`), matching both siblings.

```r
find_specific_markers_sc <- S7::new_generic(
  name = "find_specific_markers_sc",
  dispatch_args = "object",
  fun = function(
    object,
    column_of_interest,
    reference_group = NULL,
    method = "wilcox",
    alternative = c("greater", "less", "twosided"),
    min_prop = 0.05,
    downsampling = TRUE,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)
```

`reference_group = NULL` loops over every level of the column; naming a group
runs that one arm only.

Method body:

1. `match.arg(alternative)`, then the same checkmate block as
   `find_all_markers_sc` (lines 191-201) plus
   `checkmate::qassert(reference_group, c("S1", "0"))`.
2. `obs_data <- object[[c("cell_id", column_of_interest)]][!is.na(get(column_of_interest))]`,
   then `unique_groups <- sort(unique(obs_data[[column_of_interest]]))`. Sorting
   matters: the downsample seed is derived from a group's position in this
   vector, so the same group gets the same subsample in every arm.
3. Hard error via checkmate if fewer than 2 groups. Warn (do not error) at
   exactly 2: one-vs-many with a single rival is the pairwise case, and upstream
   asserts the two agree exactly.
4. `reference_group` not `NULL` -> `checkmate::assertChoice(reference_group, unique_groups)`.
5. Build the barcode list per group once, `purrr::map`, and downsample each
   group independently to 100,000 cells with `set.seed(seed + group_position)`
   when `downsampling` is `TRUE`. Same cap as the sibling, applied per group
   rather than to the pooled complement.
6. Loop the reference arms. Per arm: reference barcodes from the list, rivals
   are `unique_groups` minus the reference in that same order, then

   ```r
   rs_calculate_dge_one_vs_many(
     f_path = get_rust_count_cell_f_path(object),
     cell_indices_ref = get_cell_indices(object, cells_ref, rust_index = TRUE),
     cell_indices_other = purrr::map(
       rival_cells,
       \(x) get_cell_indices(object, x, rust_index = TRUE)
     ),
     min_prop = min_prop,
     alternative = alternative,
     verbose = 0L
   )
   ```

   `rust_index = TRUE` gives the 0-based indices the binding expects.
7. Gene ids via `get_gene_names(object)[res$genes_to_keep]`, recycled
   `n_rivals` times for the per-comparison table (rival-major, matching the
   binding's flattening order). Progress `message()` per arm gated on
   `.verbose`, as in `find_all_markers_sc`.
8. Skip the arm with a warning if `sum(res$genes_to_keep) == 0` rather than
   binding zero-row tables silently.
9. `data.table::rbindlist` both tables across arms, `setcolorder`, hand to
   `new_sc_specific_markers()`.

Roxygen: full block on the generic with `@param` for every argument, `@return`
describing both tables via `\itemize{}`, `@export`, and
`@references Lun, et al., F1000Research, 2016` plus
`Soneson and Robinson, Nat Methods, 2018`. The method gets
`#' @method find_specific_markers_sc ScOrScSubset` only, matching the siblings.

## Step 6: tests

`inst/tinytest/test_sc_processing.R`, new `### find specific markers` block
after the find-all-markers block (ends line 1086). Reuse the existing fixture:
`generate_single_cell_test_data()` gives 3 cell types in `cell_grp` with
disjoint marker blocks (`cell_type_1` -> `gene_001`..`gene_010`, `cell_type_2`
-> `gene_011`..`gene_020`, `cell_type_3` -> `gene_021`..`gene_030`). That is
exactly the specificity structure this function is meant to resolve.

Cover:

- Class is `ScSpecificMarkers`; both slots are `data.table`s with the expected
  columns (`checkmate::testDataTable`).
- Looping form: `unique(summary$ref_grp)` is the 3 cell types;
  `nrow(per_comparison) == 2 * nrow(summary)`.
- Specificity, the point of the whole thing: for each reference, its own marker
  block clears `min_auroc > 0.9` and `simes_fdr <= 0.05`, and the other two
  blocks do not clear `min_auroc > 0.9`. This is what a pooled one-vs-rest
  cannot guarantee.
- `worst_rival` is always one of the other two group names, never the reference.
- `min_rank >= 1L`.
- Single-reference form: only that `ref_grp` comes back, and its summary rows
  are identical to the looping form's rows for that group (same seed, same
  subsample, so exact equality holds).
- Agreement with the sibling: on a two-level column, `per_comparison$auroc`
  equals `find_markers_sc`'s `auroc` for the same two groups. Upstream's
  `test_one_vs_many_matches_pairwise` asserts this, and the gene filter is the
  same set when there is only one rival.
- Errors: `column_of_interest` not in obs, `reference_group` not a level of the
  column. Check the actual condition message rather than a bare
  `expect_error()`, which passes for the wrong reason.

## Step 7: docs and registration

- `devtools::document()` regenerates `R/extendr-wrappers.R` (never hand-edit),
  `NAMESPACE` and `man/`.
- `_pkgdown.yml`: add `find_specific_markers_sc` after `find_all_markers_sc`
  (line 521) in "Single cell analysis methods". Do **not** index
  `rs_calculate_dge_one_vs_many`: `rs_calculate_dge_mann_whitney` is not in
  pkgdown either, since the on-disk streaming `rs_*` functions are left out.
- `R/zzz.R`: `"ScSpecificMarkers"` into the print-registration vector.

## Files touched

| File | Change |
| --- | --- |
| `src/rust/Cargo.toml` | temporary path patch, revert before commit |
| `src/rust/src/single_cell/r_sc_analysis.rs` | binding + module registration |
| `R/methods_sc_dge_diffabund.R` | generic + `ScOrScSubset` method |
| `R/classes_single_cell_others.R` | `new_sc_specific_markers` + print |
| `R/zzz.R` | print registration |
| `_pkgdown.yml` | index the new generic |
| `inst/tinytest/test_sc_processing.R` | new test block |
| `R/extendr-wrappers.R`, `NAMESPACE`, `man/` | generated |

## Verification

```sh
export DEV_BUILD=true
cargo clippy --manifest-path=src/rust/Cargo.toml
air format .
```

```r
Sys.setenv(DEV_BUILD = "true")
devtools::document()
devtools::install()
tinytest::run_test_file("inst/tinytest/test_sc_processing.R")
```

Then a manual smoke check that the summaries actually separate the marker
blocks, not just that the shapes line up:

```r
res <- find_specific_markers_sc(sc_object, "cell_grp", .verbose = FALSE)
res
res$summary[ref_grp == "cell_type_1"][order(-median_auroc)][1:15]
res$per_comparison[ref_grp == "cell_type_1" & gene_id == "gene_001"]
```

Genes 1-10 should top `cell_type_1` on `median_auroc` with `min_auroc` near 1,
and `worst_rival` should point at whichever of the other two types shares the
most signal.

Finally revert `src/rust/Cargo.toml` to the git/version patch before committing.
