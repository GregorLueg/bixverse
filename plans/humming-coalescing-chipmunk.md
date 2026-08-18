# Consensus NMF: crate gaps and R wiring

## Context

`bixverse-rs` on branch `feat-nmf-clustering` (v0.4.5, unpublished) gained a
cNMF implementation (Kotliar et al. style): pool the components of `n_runs`
HALS restarts, drop low-density outliers, k-means the survivors, take the
per-cluster median as the consensus factor, refit the partner factor against
it. Plus a k-sweep that reports stability against reconstruction error so you
can pick `k`.

Six entry points exist in the crate: `nmf_consensus_run` / `nmf_k_sweep_run`
(dense), `*_sparse`, and `*_sc` (disk-backed). **The meta-cell variants were
never added.** The R package has none of it wired up.

Two jobs: fill the meta-cell gap in the crate, then wire consensus NMF and the
k-sweep into the R package for bulk, single cell and meta cells, with the
usual `rs_` bindings, `params_*()` wrapper, checkmate pair, S7 generics and
methods, tests and vignettes.

The single-cell vignette currently punts on this explicitly
(`vignettes/bag_of_genes_single_cells.qmd:834`, "the classical move is to
cluster these..."), so there is a ready-made hole to fill.

## Decisions

Settled with you: local path override for the crate; k-sweep returns a
detached S3 result on all three paths; expose both consensus targets with
`"h"` as the default; add a ggplot2 plot method for the sweep.

Three things I am deciding as I go, flagging because they are not obvious:

1. **Drop `clusters.consensus` at the binding.** Verified in
   `consensus.rs:791-806`: for `HRows` the consensus matrix *is* `h`
   (a `.clone()`), for `WColumns` it is `t(w)`. Returning it ships a second
   copy of a matrix the caller already has. Nothing is lost.
2. **No `consensus_seed` knob.** The crate reads a `consensus_seed` key for
   the k-means. Exposing it alongside the restart `seed` gives one call three
   seeds, and a user setting `seed = 1` would silently leave k-means at 42.
   The R methods inject `consensus_seed = seed` into the list before the call.
   One knob.
3. **`rel_error` / `rel_run_errors`, not `error` / `losses`.**
   `ConsensusNmfResult.error` and `run_errors` are divided by `||V||_F^2`
   (`consensus.rs:809-814`). `StabilisedNmfResult.losses`, already exposed via
   `get_nmf_stability()`, is absolute. Two same-shaped fields differing by
   orders of magnitude in one package is a bug report waiting to happen.

One concern on your `consensus_target` call, doing it as asked anyway: on the
single-cell path `"w"` pools a `(k * n_runs) x n_cells` dense matrix and runs
an exhaustive cosine kNN over rows of dimension `n_cells`. At k=20, n_runs=50,
200k cells that is ~800 MB before the kNN. I will expose it, document the cost
in `@details`, and emit a `warning()` from the sc/mc methods when `"w"` is
combined with a large cell count.

## Stage 0: crate, meta-cell entry points

`~/repos/shared/bixverse-rs/src/single_cell/mc_analysis/nmf_mc.rs`

Add `nmf_consensus_run_mc` and `nmf_k_sweep_run_mc` mirroring the two existing
shims line for line: CSR to CSC transpose guard under
`verbosity.detailed_verbosity()`, `Instant` timing, delegate to
`nmf_consensus_run_sparse` / `nmf_k_sweep_run_sparse`. `data:
CompressedSparseData2<f32>` by value, `consensus_params:
Option<ConsensusParams<f32>>` after `nmf_hals_params`, then `n_runs`,
`base_seed: usize`, `verbose: usize`.

The CSR guard is load-bearing, not cosmetic: `nmf_process_sparse`
(`nmf_preprocessing.rs:138-146`) indexes `indptr` per column, and
`mc_counts_to_list()` always emits `cs_type = "csr"`. Skip the transpose with
`preprocessing != "none"` and you get silent garbage, not an error.

`nmf_hals/mod.rs` declares `pub mod consensus;` with no glob re-export, so
`use crate::methods::nmf_hals::*` will not pull the types in. Import
explicitly:

```rust
use crate::methods::nmf_hals::consensus::{
    ConsensusNmfResult, ConsensusParams, KSweepEntry,
};
```

Inline `#[cfg(test)] mod tests`: CSR input with `preprocessing = "sd"` matching
the CSC equivalent (the test that actually proves the guard), `w`/`h` shapes,
and a k-sweep returning one entry per `k`.

`cargo test` and `cargo clippy --all-targets --features single-cell,multi-modal`.

## Stage 1: boundary, params, bulk binding

`src/rust/Cargo.toml`: bump `bixverse-rs` to `0.4.5` and add

```toml
[patch.crates-io]
bixverse-rs = { path = "../../../../../../bixverse-rs" }
```

That path is worktree-relative. It becomes `../../../bixverse-rs` from the main
checkout. `src/rust/Cargo.lock` is tracked, and the patch rewrites the
`bixverse-rs` entry from a registry source with a checksum to a path source
without one, which breaks CI. See Stage 6.

`R/param_wrappers.R`, after `params_nmf_hals()` (~line 472):

```r
params_nmf_consensus(
  consensus_target = c("h", "w"),
  n_neighbours = 0L,        # 0 = auto, ceil(0.3 * n_runs)
  density_threshold = 0.5,  # >= 2 disables the filter
  kmeans_iters = 100L,
  kmeans_n_init = 3L
)
```

`R/checkmate_extensions.R`, after `assertNmfHals` (~line 859):
`checkNmfConsensus` / `assertNmfConsensus` on the existing
`check_list_shape()` + `apply_qtest_rules()` + `apply_choice_rules()`
skeleton. Rules: `n_neighbours = "I1[0,)"`, `density_threshold = "N1[0,2]"`,
`kmeans_iters = "I1[1,)"`, `kmeans_n_init = "I1[1,)"`, choice
`consensus_target = c("h", "w")`.

`src/rust/src/methods/r_nmf_bulk.rs`: `rs_nmf_consensus_bulk` and
`rs_nmf_k_sweep_bulk`, added to the file-local `extendr_module!` (no `lib.rs`
change). `@export` plus `@keywords internal`, `@description` opening with the
lifecycle badge, mirroring `rs_nmf_single_bulk`.

Flat return, R-side constructors do the reshaping as everywhere else:

```
w, h, rel_error, rel_run_errors, labels, local_density, kept,
silhouette, stability, cluster_sizes, n_dropped, n_empty_clusters
```

`labels` is `Vec<Option<usize>>` mapping to integer with `NA` for dropped;
`labels` and `kept` shifted to 1-based; `f32` scalars cast to `f64`.

k-sweep returns columnar vectors that R turns straight into a data.table:

```
k, stability, best_error, median_error, consensus_failed,
n_dropped, n_empty_clusters, n_converged
```

`rextendr::document()` with `DEV_BUILD=true`, then a scratch test calling
`rs_nmf_consensus_bulk` directly on a small matrix: `length(labels) == k *
n_runs`, `sum(!is.na(labels)) == length(kept)`, `length(silhouette) ==
length(kept)`. Debug the params plumbing on the fast dense f64 path before any
class machinery exists.

## Stage 2: bulk user API

`R/helpers.R`: `.nmf_cluster_dt()` internal, building the per-component
diagnostics data.table from the flat vectors. Shared by all three paths. One
row per pooled component, with `component_id` as `run_%02d.comp_%02d`.

That naming is not arbitrary: `pool_components()` (`consensus.rs:279`) fills
pooled row `r * k + j` from restart `r` component `j`, which is exactly the
order `new_stabilised_nmf_result()` already builds `w_all` colnames in. So the
diagnostics join straight onto `get_w()` of a stabilised fit.

Columns: `component_id`, `run`, `component`, `pooled_idx`, `cluster` (`NA` if
dropped), `local_density` (`Inf` if collapsed), `silhouette` (`NA` if
dropped), `kept` (logical). One joinable object rather than five ragged
vectors of three different lengths.

`R/methods_coexp_nmf.R`:

- `consensus_nmf_bulk` generic + `BulkCoExp` method, following `nmf_bulk` line
  for line: `match.arg`, checkmate block plus `assertNmfConsensus`,
  `.assert_bulk_detection_method(..., allow_unset = TRUE)` early return,
  `.get_bulk_target_mat()`, negativity check, `consensus_seed` injection,
  `rs_` call, `sample_activity <- res$w` / `gene_loadings <- t(res$h)`,
  `sprintf("comp_%02i", seq_len(k))` dimnames, `.nmf_modules_from_w()`, then
  `new_bulk_module_result(method = "nmf-based", ...)` with the consensus
  material in `diagnostics` (`stability`, `rel_error`, `rel_run_errors`,
  `n_dropped`, `n_empty_clusters`, `clusters`). `params$nmf_fit` gets
  `consensus = TRUE`, `n_runs`, `nmf_consensus_params`.

  Reuse `"nmf-based"` rather than adding a method string:
  `new_bulk_module_result()`'s `method` is a closed `assertChoice` set
  (`R/classes_bulk.R:1166`) and every getter keys on
  `inherits(final_results, "BulkModuleResult")`.

- `nmf_k_sweep_bulk` generic + method, returning a detached
  `NmfKSweepResult`. Does not set `detection_method` and does not touch
  `final_results`: locking a purely diagnostic call out of every other module
  method would be maddening.

- Extend `get_nmf_stability()`. It currently early-returns `NULL` with a
  warning when `diagnostics[["losses"]]` is absent, and a consensus fit has
  `rel_run_errors`, not `losses`. Without this it silently claims a fit that
  exists does not.

Pre-validate in R before minutes of restarts: `qassert(k, "I1[2,)")`,
`qassert(n_runs, "I1[2,)")`. The crate enforces both but only after the work.

Error surface: let `NmfConsensusTooFewComponents` and
`NmfConsensusEmptyCluster` propagate, they are data-dependent and only knowable
after the restarts. Wrap the `rs_` call in `tryCatch` and re-`stop()` with the
original message plus one sentence naming the escape hatches
(`params_nmf_consensus(density_threshold = 2)`, raise `n_runs`, run the sweep
first). No catch-and-retry and no degraded partial result: an empty cluster
means a zero row in the consensus factor, and instability is the signal cNMF
is built to report.

Tests in `inst/tinytest/test_bulk_nmf.R`, same
`synthetic_bulk_cor_matrix(generator = "non_negative_factor")` fixture: recall
against planted modules at `k = 5`, `n_runs = 10`; `stability` within
`[-1, 1]`; sweep over `k = 2:5` giving 4 ordered rows; same-seed
reproducibility; `n_runs = 1` errors; a deliberately tight
`density_threshold` erroring with a message naming the threshold.

Use `density_threshold = 2` (filter off) as the test default.
`resolve_n_neighbours` is `ceil(0.3 * n_runs)`, so at small `n_runs` a 0.5
threshold on a noisy fixture can leave fewer than `k` survivors, which is a
hard error, which is a flaky test.

## Stage 3: single cell

`src/rust/src/single_cell/r_sc_analysis.rs` after `rs_nmf_multi_sc` (~line
1745): `rs_nmf_consensus_sc`, `rs_nmf_k_sweep_sc`. Same `f_path_gene` /
`&[i32]` / `.r_int_convert()` / `ParallelSparseReader::new` pattern, f32,
`as f64` on scalars, `.r_float_convert()` on vectors.

Validate `k_range` R-side: `nmf_k_sweep_run_sc` loads the matrix from disk
before `nmf_k_sweep_run_sparse` validates, so a typo costs a full load.

`R/classes_single_cell_others.R` after `StabilisedNmfResult` (~line 3300):

- `new_consensus_nmf_result(nmf_res, gene_ids, cell_ids, cell_indices,
  source_class, params)` producing `ConsensusNmfResult`. Applies the same flip
  as `new_nmf_result()`: `w <- t(nmf_res$h)` (genes x k), `h <- t(nmf_res$w)`
  (k x cells). On the sc/mc path V is cells x genes, so Rust `w` is cells x k
  and `h` is k x genes, despite what the existing binding doc comment on
  `rs_nmf_single_sc` claims. Fix that doc comment while in there.
  Carries `$clusters` (the data.table), `$stability`, `$rel_error`,
  `$rel_run_errors`, `$cluster_sizes`, `$n_dropped`, `$n_empty_clusters`.
- `print.ConsensusNmfResult`, `dim`, `as.matrix`, `get_w`, `get_h`,
  `get_params`, `get_data` (mirroring `get_data.NmfResult` so consensus
  activations flow into `add_sc_new_obs()`), and a `get_stability()` generic.
- `new_nmf_k_sweep_result()` returning a **subclassed data.table**
  (`c("NmfKSweepResult", "data.table", "data.frame")`) with params in
  `attr()`. `[`, `order`, `melt` and ggplot keep working and
  `print.NmfKSweepResult` is a header plus `NextMethod()`. Much less code than
  getters over a wrapper.
- `plot.NmfKSweepResult`: stability and relative error against `k`, ggplot2
  (already an Import). Exported like `plot.synthetic_matrix_simple`, no zzz
  entry needed. Covers all three paths since the class is shared.

Build the stability column as `data.table::fifelse(consensus_failed,
NA_real_, stability)` rather than trusting `f32::NAN` to survive the boundary
as something `is.na()` likes, and `warning()` with the count of failed `k`
values so nobody plots a curve with silent holes. Warn separately when any row
has `n_empty_clusters > 0`: `silhouette_cosine_unit` returns 0 when fewer than
two clusters are populated (`clustering_metrics.rs:139-142`), so a stability of
0 in the sweep is ambiguous.

`R/zzz.R`: add `"ConsensusNmfResult"` and `"NmfKSweepResult"` to the `classes`
vector for manual print registration.

`R/base_generics_sc.R` after `stabilised_nmf_sc` (~line 1205):
`consensus_nmf_sc` and `nmf_k_sweep_sc` generics, `@inheritParams nmf_sc`.

`R/methods_sc_pathways.R` after the `stabilised_nmf_sc` method (~line 1085):
both `SingleCells` methods, reusing `.resolve_sc_nmf_selection()` unchanged.

Put a rough memory formula in `@details`. `nmf_consensus` holds the
`SparseInput` (CSR plus CSC copies), `w_all` (`n_cells x k*n_runs` dense) and
`h_per_run` (`n_runs x k x n_genes` dense) at once. Already true of
`stabilised_nmf_sc`, but consensus is what makes people want `n_runs = 100`.
Steer toward the `MetaCells` path.

Tests in `inst/tinytest/test_sc_gs_activity.R` after the existing NMF block,
using `sc_test_object()` / `sc_test_prepped()`, `k = 3`, `n_runs = 5`,
`density_threshold = 2`: `W` is `(n_hvg, k)`, `H` is `(k, n_cells_kept)`, both
named and non-negative, `nrow(clusters) == k * n_runs`, `component_id` matches
`colnames(get_w(stabilised_fit))`, same-seed reproducibility, `n_runs = 1`
errors, sweep over `k = 2:4` gives 3 rows.

## Stage 4: meta cells

`src/rust/src/meta_cell/r_mc_analysis.rs` after `rs_nmf_multi_mc` (~line 545):
`rs_nmf_consensus_mc`, `rs_nmf_k_sweep_mc`, via `list_to_sparse_matrix(...,
true)` then `cast_compressed_sparse_data_f32`.

Add `@keywords internal` to these two and to the existing `rs_nmf_single_mc` /
`rs_nmf_multi_mc`, which lack it and are therefore the only NMF bindings listed
in `_pkgdown.yml`. Then delete `_pkgdown.yml:777-778`. Both halves have to move
together: pkgdown fails on exported-but-unindexed topics unless
`@keywords internal` is present.

`R/methods_mc_analysis.R` after the `stabilised_nmf_sc` `MetaCells` method
(~line 560): both methods, reusing `.resolve_mc_nmf_selection()` and
`mc_counts_to_list()` unchanged, passing `cell_indices = sel$cell_indices_rust`
into the constructor.

Tests in `inst/tinytest/test_mc_analysis.R` after the existing NMF block
(~line 640), mirroring the sc ones on a `MetaCells` object.

## Stage 5: docs

`_pkgdown.yml`: `consensus_nmf_bulk`, `nmf_k_sweep_bulk` into the bulk
co-expression block (~197-231); `consensus_nmf_sc`, `nmf_k_sweep_sc`,
`get_stability`, `params_nmf_consensus` into the single cell block (~548-566);
delete lines 777-778.

`NEWS.md` under a new heading, `DESCRIPTION` version bump, `src/rust/Cargo.toml`
version kept in sync.

Vignettes last, once the API is frozen:

- `vignettes/bulk_coexpression_modules.qmd` section "Method 5" (420-494): a
  `k = 2:6` sweep with the plot, then `consensus_nmf_bulk()` at the chosen
  `k`, folded into the existing planted-module recovery comparison.
- `vignettes/bag_of_genes_single_cells.qmd` (674-841): replace the
  "Consensus across runs" bullet that currently punts with a real subsection.
- `vignettes/meta_cells.qmd` section at line 672: short, framed as the path
  where consensus NMF is actually affordable.

Vignettes build on a 2-core, 8 GB runner. Keep `n_runs <= 10` and `k` small,
and say in prose that these numbers are too small for real work.

## Stage 6: un-patch before merge

The path override cannot be merged. Before the PR lands: publish `bixverse-rs`
(0.4.5, or 0.4.6 if the meta-cell entry points go out as a separate release),
drop the `[patch.crates-io]` block, point at the published version, then
`cargo update --manifest-path=src/rust/Cargo.toml -p bixverse-rs` and confirm
`src/rust/Cargo.lock` shows a `registry+https://...` source with a checksum
again.

I will flag this rather than commit the patch block silently. Say the word if
you would rather I keep the patch out of git entirely and you apply it locally.

## Verification

```sh
# crate
cargo test --manifest-path=$HOME/repos/shared/bixverse-rs/Cargo.toml \
  --features single-cell,multi-modal
cargo clippy --manifest-path=$HOME/repos/shared/bixverse-rs/Cargo.toml \
  --all-targets --features single-cell,multi-modal

# R package
cargo clippy --manifest-path=src/rust/Cargo.toml
air format .
```

```r
Sys.setenv(DEV_BUILD = "true")
devtools::document()
devtools::install()

tinytest::run_test_file("inst/tinytest/test_bulk_nmf.R")
tinytest::run_test_file("inst/tinytest/test_sc_gs_activity.R")
tinytest::run_test_file("inst/tinytest/test_mc_analysis.R")
tinytest::test_package("bixverse")

quarto::quarto_render("vignettes/bulk_coexpression_modules.qmd")
```

End-to-end sanity on the bulk path, which is fast enough to iterate on:

```r
synth <- synthetic_bulk_cor_matrix(
  params_synthetic_bulk_rnaseq(
    num_samples = 60L, num_genes = 300L,
    module_sizes = c(50L, 50L, 50L),
    generator = "non_negative_factor", seed = 123L
  )
)
obj <- BulkCoExp(
  raw_data = t(synth$counts),
  meta_data = data.table::data.table(sample_id = colnames(synth$counts))
) |>
  preprocess_bulk_coexp(scaling = FALSE, hvg = NULL, .verbose = FALSE)

sweep <- nmf_k_sweep_bulk(obj, k_range = 2:6, n_runs = 10L, .verbose = FALSE)
plot(sweep)

fit <- consensus_nmf_bulk(obj, k = 3L, n_runs = 10L, .verbose = FALSE)
get_nmf_stability(fit)
get_nmf_modules(fit)[, .N, by = module_id]
```

Recall against `synth$module_data` should land in the same ballpark as
`nmf_bulk()` (0.7 to 0.83 on this generator) and ideally above it, since that
is the point of the consensus step.

## Risks

- **CI breaks on the path patch.** Handled by Stage 6, but it has to actually
  happen before merge.
- **Flaky consensus tests.** The density filter erroring on small fixtures is
  the main hazard. Mitigated by `density_threshold = 2` in tests, with one
  deliberate test for the error path.
- **Vignette build time.** Consensus is `n_runs` full NMF fits, and the sweep
  is that times `length(k_range)`. The 2-core runner is the binding
  constraint.
