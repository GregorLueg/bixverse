# Palantir + PAGA: the R side

## Context

Palantir and PAGA landed in `bixverse-rs` (local checkout,
`/Users/gregorlueg/repos/shared/bixverse-rs/src/single_cell/sc_trajectory/`,
crate version 0.4.1). Nothing in the R package can reach them yet. This adds the
full R surface: extendr bindings, a parameter wrapper, checkmate extensions,
two S3 result classes and S7 generics that dispatch across `SingleCells`,
`SingleCellsSubset`, `MetaCells` and `SingleCellsMultiModal`.

Both algorithms consume a kNN graph, so every single-cell class in the package
can feed them. The heavy lifting is already done upstream, including
`PalantirParams::from_r_list` (`bixverse-rs/src/single_cell/sc_r_wrappers.rs:3132`),
so this repo only needs the thin binding plus wrappers.

## Decisions taken

- Results come back as new S3 classes `PalantirRes` and `PagaRes`, returned
  rather than written into the object. Matches `SingleCellFastClusters` and
  `ScMiloRRes`.
- Both methods take `modality = c("rna", "adt", "wnn")`. PAGA reads only kNN
  indices so it is safe on any graph. Palantir on WNN warns once, because
  `rs_wnn` stamps `dist_metric = "kernelised pseudo-distance"`, which is not a
  metric and makes `squared_dist = FALSE`.
- `early_cell` and `terminal_states` are cell names, resolved against
  `get_knn_obj(object)$used_cells`. kNN indices are positions in `used_cells`,
  not global cell indices, so this is also how output indices map back.

## What already exists and gets reused

- `knn_data_to_rust()` (`src/rust/src/single_cell/utils.rs:132`) turns a
  `SingleCellNearestNeighbour` list into `(Vec<Vec<usize>>, Vec<Vec<f32>>, k,
  dist_metric)`, exactly what `run_palantir` wants. `knn_indices_processing()`
  (same file, line 79) covers PAGA.
- `rs_metacell_density` (`src/rust/src/single_cell/r_sc_metacells.rs:935`) is a
  line-for-line template: same kNN unpack, same `squared_dist = distance ==
  "euclidean"`, same `KnnParams::from_r_list`, same verbosity plumbing.
- `sparse_data_to_list` (Rust) pairs with `sparse_list_to_mat`
  (`R/helpers.R:270`) for the PAGA connectivity matrices. Both agree on the
  `cs_type` key.
- `faer_to_r_matrix` handles `Mat<f32>` directly via `FaerRType`. No manual
  conversion needed for `branch_probs` / `multiscale`.
- `checkKnnParams` + `KNN_PARAM_NAMES` (`R/checkmate_extensions.R:7,1308`),
  `apply_qtest_rules`, `apply_choice_rules`, `check_list_shape`.
- `params_knn_defaults()` (`R/param_constructors.R:46`), `parse_verbosity()`
  (`R/functions_single_cell.R:24`).

## Step 1: point at the local crate

`src/rust/Cargo.toml` currently patches `bixverse-rs` to a git branch. Swap the
`[patch.crates-io]` entry to
`bixverse-rs = { path = "/Users/gregorlueg/repos/shared/bixverse-rs" }`. The
declared dep stays `version = "0.4.0"`; local 0.4.1 satisfies it. Leave the
`bixverse` crate version at 0.4.7 until we are ready to cut a release.

## Step 2: Rust bindings

New file `src/rust/src/single_cell/r_sc_trajectory.rs`, module doc with `//!`,
`extendr_module!` block, registered in `src/rust/src/single_cell/mod.rs`, and in
`src/rust/src/lib.rs` (both the `pub use` list and the `extendr_module!` block).

### `rs_palantir`

```rust
fn rs_palantir(
    knn_data: List,
    palantir_params: List,
    early_cell: usize,          // 0-indexed position in used_cells
    terminal_states: Nullable<Vec<i32>>,  // 0-indexed, or NULL
    seed: usize,
    verbose: usize,
) -> Result<List>
```

Body: `knn_data_to_rust` -> `squared_dist = dist_metric == "euclidean"` ->
`PalantirParams::from_r_list` -> `run_palantir(...)` -> `.to_extendr()?`.

Returns a flat list: `pseudotime`, `entropy` (both `r_float_convert()` on the
`Vec<f32>`), `branch_probs` and `multiscale` via `faer_to_r_matrix`,
`terminal_states`, `waypoints`, `start_cell` as `i32` (still 0-indexed, R side
shifts), plus `iterations`, `converged`, `repair_edges`, `stranded_waypoints`.

### `rs_paga`

```rust
fn rs_paga(
    knn_mat: RMatrix<i32>,
    partitions: Vec<i32>,   // 0-indexed cluster labels
    n_partitions: usize,
) -> Result<List>
```

Takes the raw kNN matrix rather than the full kNN list, since PAGA ignores
distances. `knn_indices_processing` -> `run_paga::<f64>` -> returns
`connectivities`, `connectivities_tree` (both `sparse_data_to_list`) and
`sizes`.

Roxygen on both in the extendr doc-comment style used across
`r_sc_metacells.rs`, with `lifecycle::badge("experimental")` and
`@keywords internal`. References: Setty et al., Nat. Biotechnol., 2019;
Wolf et al., Genome Biol., 2019.

## Step 3: parameter wrapper and checkmate extension

`params_sc_palantir()` in `R/param_constructors.R`, flat list mirroring
`PalantirParams::from_r_list` exactly. Follows the `params_sc_neighbours()`
shape: named args, then `modifyList(params_knn_defaults(), knn, keep.null =
TRUE)` folded in with `c()`.

Fields and defaults, matching the Rust `Default`: `n_dcs = 10L`,
`n_eigs = NULL`, `knn = 30L`, `num_waypoints = 1200L`,
`scale_components = TRUE`, `use_early_cell_as_start = FALSE`,
`max_iterations = 25L`, `branch_prob_threshold = 0.01`, plus the advanced
Lanczos knobs `lanczos_basis_size = NULL`, `lanczos_max_restarts`,
`lanczos_tol`, and the kNN block.

Naming footgun to document: `knn` (self-inclusive neighbour count for the
multiscale search, Palantir's own parameter) is distinct from `k` (the kNN
params block). They do not collide because the Rust side reads them under
different keys, but the roxygen must say so. Also note `ann_dist` and `k` from
the kNN block are overridden inside `geodesic_graph`.

Paired `checkScPalantirParams` / `assertScPalantirParams` in
`R/checkmate_extensions.R`, in the single-cell section. Structure copies
`checkScNeighbours` (line 2272): `check_list_shape`, then
`checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])`, then `apply_qtest_rules`
for the Palantir-specific fields. Range rules worth enforcing, taken from the
Rust constants: `n_dcs = "I1[3,)"`, `knn = "I1[6,)"` (`MIN_KNN`),
`max_iterations = "I1[2,)"` (`MIN_MAX_ITERATIONS`),
`branch_prob_threshold = "N1[0,1]"`, `num_waypoints = "I1[1,)"`,
`n_eigs = c("0", "I1[3,)")`.

PAGA gets no wrapper. `run_paga` has no params struct upstream.

## Step 4: S3 result classes

Appended to `R/classes_single_cell_others.R`, which already holds the small
single-cell S3 side-cars. Two constructors plus print methods, both class names
added to the `classes` vector in `R/zzz.R`.

`new_palantir_res()` -> `PalantirRes`:

- `pseudotime` - data.table `cell_id`, `pseudotime`, `entropy`, ordered as
  `used_cells`.
- `branch_probs` - numeric matrix, rownames `used_cells`, colnames the terminal
  cell names.
- `terminal_states`, `waypoints`, `start_cell` - character cell names, mapped
  back through `used_cells` (add 1 to the Rust indices first).
- `multiscale` - numeric matrix, rownames `used_cells`.
- `run_info` - list with `iterations`, `converged`, `repair_edges`,
  `stranded_waypoints`, `modality`, `n_waypoints`.

`new_paga_res()` -> `PagaRes`:

- `connectivities`, `connectivities_tree` - via `sparse_list_to_mat`, dimnames
  set to the cluster levels.
- `sizes` - data.table `cluster`, `n_cells`.
- `params` - list with `cluster_col`, `modality`.

Print methods follow `print.SingleCellNearestNeighbour`
(`R/classes_single_cell_others.R:75`): `cat` with `sprintf`, `invisible(x)`.

## Step 5: S7 generics and methods

New file `R/methods_sc_trajectory.R`. Generics defined next to their methods,
as `fast_cluster_sc` does (`R/methods_sc_processing.R:1565`). Both dispatch on
`ScOrMc`, which covers `SingleCells`, `SingleCellsSubset`, `MetaCells` and, via
`parent = SingleCells`, `SingleCellsMultiModal`.

```r
run_palantir_sc(object, early_cell, terminal_states = NULL,
                modality = c("rna", "adt", "wnn"),
                palantir_params = params_sc_palantir(),
                seed = 42L, .verbose = TRUE)

run_paga_sc(object, cluster_col, modality = c("rna", "adt", "wnn"),
            .verbose = TRUE)
```

Method body, shared shape:

1. `match.arg(modality)`; `checkmate::assertTRUE(S7::S7_inherits(...))` over the
   three base classes; `assertScPalantirParams`; `qassert` the rest.
2. Guard `modality != "rna" && !S7_inherits(object, SingleCellsMultiModal)` with
   the same `stop()` message `find_clusters_sc` uses
   (`R/methods_sc_processing.R:1489`).
3. `knn <- get_knn_obj(object, modality = modality)`; early return with a
   warning naming `find_neighbours_sc()` if `NULL`, matching the existing
   idiom.
4. Palantir only: resolve `early_cell` and `terminal_states` against
   `knn$used_cells` via `match()`, `stop()` on any unmatched name, subtract 1 for
   Rust. Warn once when `knn$dist_metric` is neither `"euclidean"` nor
   `"cosine"` (the WNN case).
5. PAGA only: `labels <- factor(unlist(object[[cluster_col]], use.names =
   FALSE))`; assert `length(labels) == knn$n` with a message pointing at the
   kNN/obs mismatch; pass `as.integer(labels) - 1L` and `nlevels(labels)`.
   Positional alignment is the established convention here, same as how
   `find_clusters_sc` writes memberships back.
6. Call the `rs_*` function with `verbose = parse_verbosity(.verbose)`, wrap
   into the S3 class.

Full roxygen on each generic: `@description`, every `@param` typed, `@returns`
with `\itemize{}` over the S3 fields, `@export`, `@references`. `@method` tag on
the `S7::method()` assignment.

## Step 6: docs and package wiring

- `devtools::document()` regenerates `man/` and `NAMESPACE`. The `rs_*`
  wrappers in `R/extendr-wrappers.R` regenerate during compile.
- `_pkgdown.yml`: add `run_palantir_sc`, `run_paga_sc`, `params_sc_palantir` to
  the "Single cell analysis methods" section (line 497). Add the two print
  methods' classes nowhere; they are not exported symbols. Per the existing
  memory note, `rs_palantir` / `rs_paga` stay out of pkgdown since they take
  on-object kNN handles rather than plain matrices.
- No `README.md` / `NEWS.md` edits.

## Step 7: tests

New `inst/tinytest/test_sc_trajectory.R`. No public reference dataset yet, so
the tests target the binding and wrapper correctness rather than numerical
agreement with the Python Palantir.

Build a synthetic Y-shaped manifold directly as an embedding matrix (trunk plus
two branches, ~900 points in ~10D with noise, `rownames` as fake barcodes), run
`generate_sc_knn()` on it with `ann_dist = "euclidean"`, then drive `rs_palantir`
and `rs_paga` directly.

Palantir assertions:

- `pseudotime` in `[0, 1]`, start cell exactly 0.
- Spearman correlation of pseudotime against the true trunk position > 0.9.
- Two terminal states detected, each within the tip region of a branch.
- `branch_probs` rows in `[0, 1]`, entropy near 0 at the tips and higher at the
  fork.
- Same seed twice gives identical output.

PAGA assertions on a three-cluster chain A-B-C:

- `connectivities` symmetric, zero diagonal, `A-B` and `B-C` both exceed `A-C`.
- `connectivities_tree` has exactly `n_clusters - 1` non-zero upper-triangle
  entries.
- `sizes` sums to the cell count; empty factor levels survive via
  `n_partitions`.

Wrapper assertions:

- `params_sc_palantir()` defaults match the Rust `Default` values field by
  field.
- `assertScPalantirParams` rejects `knn = 3L`, `max_iterations = 1L`,
  `branch_prob_threshold = 1.5`.
- `run_palantir_sc` errors on an unknown `early_cell` name, and returns cell
  names (not indices) in `terminal_states` / `start_cell`.
- `run_paga_sc` errors when the cluster column length does not match `knn$n`.
- Dispatch: a `SingleCells` and a `MetaCells` object carrying the same kNN
  produce identical results.

## Verification

```sh
cargo clippy --manifest-path=src/rust/Cargo.toml
R CMD INSTALL --preclean .
air format .
```

```r
devtools::document()
tinytest::run_test_file("inst/tinytest/test_sc_trajectory.R")
tinytest::test_package("bixverse")   # nothing else should regress
```

Manual smoke test on a real object: run `find_neighbours_sc()` then
`find_clusters_sc()`, feed the result into `run_paga_sc()` and confirm the
abstracted graph reflects the cluster structure; run `run_palantir_sc()` from a
plausible start cell and check the pseudotime gradient against the UMAP.

## Deliberately out of scope

- Plotting helpers for the PAGA graph or pseudotime overlays.
- A standalone binding for `multiscale_components` (the diffusion map already
  comes back inside `PalantirResult`).
- Vignette and Zenodo dataset. Future us.
