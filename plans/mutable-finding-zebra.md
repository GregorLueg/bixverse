# Provenance stamps for `ScCache` artefacts

## Context

`ScCache` (`R/classes_single_cell_.R:46`) holds every derived artefact of the
single cell suite in memory: PCA factors/loadings/singular values, arbitrary
named embeddings (`umap`, `harmony`, `mnn`, `cca`, `rpca`, `symphony`, ...), the
kNN object and the sNN `igraph`. The cell set those artefacts were computed on
lives somewhere else entirely, in `ScMap$cells_to_keep_idx`. Nothing links the
two.

The consequence is the foot gun the package already documents in
`vignettes/thinking_single_cell.qmd:165`. `set_cells_to_keep()`
(`R/classes_single_cell_.R:1886`) rewrites the cell set and pushes it into
DuckDB, and leaves `sc_cache` completely untouched. Re-run PCA on the new cell
set, forget the kNN, and the object carries a PCA over 8 cells next to a kNN
over 10. Nothing complains. `get_pca_factors()` slaps the *new* cell names onto
the *old* matrix rows (`:1372`), meta cell generation hands mismatched indices
to Rust, and the failure mode is either a crash deep in the binding layer or,
worse, silently mis-aligned biology.

Today there is exactly one invalidation site in the package
(`SingleCellsSubset` builds a fresh cache, `R/classes_single_cell_subset.R:109`)
and a handful of ad-hoc read-time guards (`run_paga_sc` at
`R/methods_sc_trajectory.R:308`, `assign_sc_type` at
`R/methods_sc_annotations_reference.R:130`). `remove_knn()` /
`remove_snn_graph()` exist for all four classes but nothing in the package ever
calls them.

The fix: stamp every cached artefact with the cell state it was computed on plus
pointers to what it was derived from, and validate before use.

## Decisions taken

- **Stamps ride as an attribute on the payload itself**, `attr(x,
  "bixverse_stamp")`. `ScCache` keeps its six slots and its structure is
  completely unchanged. No parallel bookkeeping means no way for stamps and
  payloads to drift apart: `remove_knn()` deletes the payload and the stamp
  disappears with it, `set_embedding()` overwriting a name replaces both.
- **Stamp shape**: `list(cells = <hash>, id = <hash>, from = <character>)`. The
  `cells` field catches cell-set drift; the `id`/`from` chain catches "PCA
  re-run on the same cells with different HVGs or `n_pcs`, kNN never
  regenerated". `from` is an **active staleness axis**, checked recursively.
- **Two severity tiers.** The S7 getters warn. Functions that hand indices to
  Rust or write results back into the object call `assert_sc_state()` at entry
  and hard-error. That matches how everything else in the package validates its
  inputs, and it keeps the getters usable as presence probes.
- **Legacy objects**: an artefact with no stamp attribute is `unknown` and
  passes. This needs no code, it falls out of `attr()` returning `NULL`.
- **`set_cells_to_keep()`**: warns listing what it just invalidated, keeps the
  artefacts. Binning a 40 minute Harmony run silently is worse than a warning.
- **New `reset_cells_to_keep()`**: the escape hatch back to a pristine object.
  Restores every cell in the binary and wipes the cache completely. With
  `force = FALSE` (the default) it will not do that behind your back, it asks.

## Design

### The stamp

```r
list(cells = <chr1>, id = <chr1>, from = <chr>)
```

- `cells` - `rlang::hash()` of the owning object's cell state at write time.
  Order sensitive: `set_cells_to_keep.ScMap` does not sort, and PCA rows come
  back from Rust in exactly `cell_indices` order, so row order is what aligns a
  matrix to a cell list.
- `id` - 16 chars, `rlang::hash(list(Sys.time(), Sys.getpid(), <counter>))` with
  a monotonic counter in a package-local env. **Not** `runif()`/`sample()`:
  consuming RNG state would shift results in every `set.seed()`-based test under
  `inst/tinytest/`.
- `from` - zero or more ids of the artefacts this one was derived from.
  `character(0)` for a root. A vector, because the WNN kNN is fused from two
  embeddings living in two different caches.

`rlang` is already in `Imports`, so no new dependency.

### Where the attribute goes

| artefact | payload | notes |
|---|---|---|
| PCA | `pca_factors` only | Factors, loadings and singular values are produced atomically by one producer; loadings are gene-space and singular values dimensionless, so a cell hash on them is meaningless. One stamp for the triple |
| embedding | `other_embeddings[[name]]` | |
| kNN | the `SingleCellNearestNeighbour` list | classed list, attribute is invisible to its `print` method |
| sNN | the `igraph` | set once, after `igraph::add_edges()` has finished building it |

### Freshness

An artefact is stale when either:

1. `stamp$cells != .sc_state_hash(x)` - the cell set moved, or
2. any id in `stamp$from` no longer exists among the object's stamped artefacts
   - the upstream was recomputed (new `id`) or removed, or
3. the upstream artefact is itself stale - recursive, depth at most
   `pca -> embedding -> knn -> snn`.

Recursion matters: without it, re-running PCA marks the kNN stale but the sNN
still resolves its parent (the untouched kNN) and passes.

A `NULL` payload has no stamp, so absent artefacts stay `unknown` and the
existing presence probes keep working: `harmony_sc`
(`R/methods_sc_batch_corr.R:898`), `fast_mnn_sc:805`, `bbknn_sc:619`, and in
`bixverse.gpu` `R/sc_gpu.R:626` and `:747`.

### Attribute mechanics, honestly

- **Stripped on the way out.** The S7 getters `remove_attr()` before returning,
  so the stamp never reaches user code, never shows up in `str()`, and never
  makes `identical()` on two numerically equal matrices return `FALSE`. Existing
  tests never compare getter output attribute-wise (they subset, which drops
  attributes anyway, or run `checkmate::testMatrix`), but stripping removes the
  category. Internal `from` resolution therefore reads the **raw cache slot**,
  not the public getter.
- **Subsetting drops attributes.** `find_neighbours_sc:1414` does `embd <-
  embd[, 1:no_embd_to_use]`, as does `bixverse.gpu` at `R/sc_gpu.R:126, 255,
  930, 1114`. Harmless: those are local copies made *after* the read, and the
  cached object keeps its attribute.
- **igraph.** The stamp is set once on a finished graph and only ever read back,
  never passed through an igraph function that would rebuild it.
- **Serialisation.** `qs2` and `rds` both preserve attributes, so
  `save_sc_exp_to_disk()` / `load_existing()` need **no change at all**.

## Implementation

### New file `R/functions_sc_cache_state.R`

Internal unless marked. This is a much smaller surface than a parallel stamp
registry would need: no key routing, no nested `stamps$embeddings[[name]]`, no
virtual-`"pca"` mapping, no skeleton normalisation.

| function | purpose |
|---|---|
| `.new_sc_id()` | 16-char id, time + PID + counter |
| `.set_stamp(obj, stamp)` | `attr(obj, "bixverse_stamp") <- stamp`; returns `obj` |
| `.read_stamp(obj)` | the stamp or `NULL` |
| `.drop_stamp(obj)` | used by the getters before returning |
| `.sc_payload(x, artefact, name, modality)` | raw cache slot read, no checking, no stripping. The only thing `from` resolution uses |
| `.sc_stamped_artefacts(x)` | walks every cache the object owns (`sc_cache`, and for multimodal `adt_cache`, `atac_cache`, `other_data$wnn`) and returns modality-qualified label + stamp pairs. Backs `from` resolution, the id lookup and the report |
| `.sc_state_hash(x)` | S7 generic, see below |
| `.stamp_ids_from(x, names, modality)` | resolves parent artefact *names* (`"pca"`, `"umap"`, `"adt:pca"`) to ids |
| `.stamp_artefact(x, artefact, name, modality, from)` | reads the raw slot, stamps it, writes it back. Returns `x` |
| `.sc_state_report(x)` | the engine. `data.table(modality, artefact, name, present, stamped, stale, reason, id, from)`. Computes `.sc_state_hash(x)` **once** and reuses it across rows |
| `.warn_sc_state(x, artefact, name, modality)` | the soft tier used by the getters |
| `check_sc_state(x, artefacts, modality)` | **exported**. `TRUE` or a message string, checkmate convention |
| `assert_sc_state` | **exported**. `checkmate::makeAssertionFunction(check_sc_state)`. The hard tier |
| `get_sc_cache_status(object)` | **exported**. The report as a `data.table`, for eyeballing |

Naming is snake_case rather than the `checkXxx`/`assertXxx` camelCase of
`R/checkmate_extensions.R`, because this validates object state rather than a
`params_*()` list, and it sits next to `get_sc_cache_status()` and
`set_cells_to_keep()` in the user's mental model. Worth a sentence of
justification in the roxygen.

Producers never touch ids. They pass parent *names* and `.stamp_artefact()`
resolves them, which keeps the churn at each producer site to one argument.

### `.sc_state_hash()`

S7 generic declared in `R/base_generics_sc.R` next to the other cache generics
(`@noRd`, unexported); methods next to each class.

- `SingleCells` (near `classes_single_cell_.R:1236`) -
  `rlang::hash(get_cells_to_keep(x))`.
- `SingleCellsMultiModal` - **no method**, inherits. All three caches and the
  wnn slot are aligned to the same kept-cell set.
- `SingleCellsSubset` (near `classes_single_cell_subset.R:343`) - same
  expression; its `get_cells_to_keep` falls back to `subset_to_original - 1L`,
  so it is defined even when `cells_to_keep_idx` is empty.
- `MetaCells` (near `classes_meta_cell.R:532`) -
  `rlang::hash(S7::prop(x, "obs_table")[["meta_cell_id"]])`. That is the
  identity `get_pca_factors.MetaCells:800` and `get_embedding.MetaCells:869`
  already use for rownames, and `merge_meta_cells`
  (`R/functions_meta_cell.R:237`) rewrites it with source prefixes so merged
  objects hash differently. Adding an obs column via `[[<-` does not change it,
  which is correct: a new cluster column must not invalidate the PCA.

Always hash a plain atomic vector, never a `data.table` - the
`.internal.selfref` pointer makes that hash unstable.

Cost is a xxhash128 over roughly 4 MB at 1M cells, single-digit milliseconds,
once per check. `.sc_state_report()` hashes once for the whole table. Do **not**
memoise onto the object; S7 has value semantics and the memo would itself need
invalidating.

### Cache routing

`.sc_payload()` and `.sc_stamped_artefacts()` need a modality router. Test
`S7::S7_inherits(x, SingleCellsMultiModal)` **first**, it has `parent =
SingleCells`. Reuse the existing `.cache_slot_from_modality()`
(`classes_single_cell_multimodal.R:329`) rather than duplicating the switch. For
`modality == "wnn"` read `other_data[["wnn"]]` tolerantly; do **not** reuse
`.integration_slot():338`, which errors.

### Leaf `set_*.ScCache` / `get_*.ScCache` methods

**Unchanged**, and so are the nine generics in `R/base_generics_sc.R:384-476`.
The leaf has no access to the cell state, and these are exported S3 methods
anyone can call on a bare `ScCache`.

One small unrelated fix worth taking while in there: `set_embedding.ScCache:179`
should reject `name == "pca"`. `get_embedding.ScCache:450` treats `"pca"` as a
virtual name that wins over `other_embeddings$pca`, so `umap_sc(slot_name =
"pca")` today writes an embedding that can never be read back.

### S7 forwarder layer

Four setters per class gain one statement after the payload write:

```r
S7::prop(x, "sc_cache") <- set_knn(x = S7::prop(x, "sc_cache"), knn = knn)
x <- .stamp_artefact(x, artefact = "knn", from = from)
return(x)
```

`from` comes out of the `...` the forwarders already carry
(`from <- list(...)[["from"]] %||% character()`), so producers that do not pass
it get a root stamp and keep working. It must stay **after** the payload in
every signature: `bixverse.gpu` passes payloads positionally at
`R/sc_gpu.R:393, 529, 645, 674-676`.

Setter sites (`set_pca_factors`, `set_embedding`, `set_knn`, `set_snn_graph`):

- `SingleCells` `classes_single_cell_.R:1939, 2008, 2034, 2058`
- `SingleCellsSubset` `classes_single_cell_subset.R:617, ~680, 696, 711`
- `MetaCells` `classes_meta_cell.R:604, 670, 695, 717`
- `SingleCellsMultiModal` `classes_single_cell_multimodal.R:787, 826, 856, 869`
  - each already computes `modality <- match.arg(...)`; pass it through.
  `set_embedding:826` has a `"wnn"` early return, stamp inside that branch too.

`set_pca_loadings`, `set_pca_singular_vals`, `remove_knn` and `remove_snn_graph`
need **no change at all**. That is the payoff from putting the stamp on the
payload.

Six getters per class gain two lines, after the `is.null(res)` early return and
before the existing `rownames(res) <- ...`:

```r
.warn_sc_state(x, artefact = "pca")
res <- .drop_stamp(res)
```

- `SingleCells` `:1356` (`get_pca_factors`), `:1426` (`get_embedding`, map the
  virtual `"pca"`), `:1471`, `:1491`, `:1511`, `:1531`
- `SingleCellsSubset` `:750, 795, 822, 832, 842, 852`
- `MetaCells` `:779, 850, 893, 912, 931, 950`
- `SingleCellsMultiModal` `:906, 952` (both the wnn and the cache branch),
  `:1002, 1015, 1028, 1042`

`get_pca_loadings` and `get_pca_singular_val` only need `.drop_stamp()` if they
ever carry one, which they do not. `get_available_embeddings` is untouched: it
is the existence probe used by `find_neighbours_sc`, `umap_sc`, `tsne_sc`,
`phate_sc`, `calculate_batch_asw_sc` and four `bixverse.gpu` entry points.

### The hard tier - who calls `assert_sc_state()`

Rule: anything that hands indices to Rust or writes a result back into the
object asserts at entry, next to its existing checkmate calls. Anything that
only reads for plotting or a metric relies on the getter warning.

| function | file:line | asserts |
|---|---|---|
| `find_neighbours_sc` | `methods_sc_processing.R:1374` | the source embedding |
| `find_clusters_sc` | `:1468` | `"snn"` - it writes a cluster column aligned to `get_cells_to_keep()` |
| `generate_bt_meta_cells_sc` | `methods_sc_aggregations.R:72` | embedding + `"knn"` |
| `generate_seacells_sc` | `:260` | embedding + `"knn"` |
| `generate_supercells_sc` | `:392` | embedding + `"knn"` |
| `run_palantir_sc` | `methods_sc_trajectory.R:138` | `"knn"` |
| `run_paga_sc` | `:275` | `"knn"` |
| `bbknn_sc` | `methods_sc_batch_corr.R:602` | `"pca"` |
| `harmony_sc`, `harmony_v2_sc`, `fast_mnn_sc`, `seurat_cca_sc`, `seurat_rpca_sc` | `:874, 1020, 785, 1155, 1281` | `"pca"`, only on the branch that consumes it |
| `get_miloR_abundances_sc` | `methods_sc_dge_diffabund.R:363` | embedding + `"knn"` |
| `meld_sc` | `:529` | embedding + `"knn"` |
| `map_symphony_query` | `methods_sc_annotations_reference.R:487` | nothing on the query, but assert the reference's own state |
| `umap_sc`, `tsne_sc`, `phate_sc` | `methods_sc_plots_helpers.R:141, 335, 517` | the source embedding, plus `"knn"` when `use_knn` |

Read-only metrics (`calculate_kbet_sc`, `calculate_batch_asw_sc`,
`calculate_batch_lisi_sc`, the pathway scorers in `methods_sc_pathways.R`) stay
on the getter warning.

### Producer sites - what each passes as `from`

Roots, `from = character()`, passed explicitly so `git grep "from ="` finds
every producer: `calculate_pca_sc` (`methods_sc_processing.R:1329, 1358`;
`methods_sc_subset_processing.R:130, 156`; `methods_mc_processing.R:444`) and
`calculate_pca_adt_sc` (`methods_sc_multimodal.R:78`, modality `"adt"`).

| producer | file:line | `from` |
|---|---|---|
| `harmony_sc`, `harmony_v2_sc` | `methods_sc_batch_corr.R:945, 1067` | `"pca"`, read unconditionally at `:902` / `:1024` |
| `fast_mnn_sc`, `seurat_cca_sc`, `seurat_rpca_sc` | `:828, 1206, 1332` | conditional on `use_precomputed_pca` (`:805, 1183, 1309`). Hoist the condition into a local and pass `"pca"` or `character()` |
| `umap_sc`, `tsne_sc`, `phate_sc` | `methods_sc_plots_helpers.R:239, 428, 607` | `embd_to_use` resolved against `cache_modality`, plus `"knn"` when the kNN path was taken. Note the `"wnn"` asymmetry at `:193`: the parent is read from the **rna** cache while the stamp is written under `modality`. A legitimate cross-modality `from`, which is why `.sc_stamped_artefacts()` walks every cache |
| `find_neighbours_sc` | `methods_sc_processing.R:1429, 1458` | kNN gets `embd_to_use`; sNN gets `"knn"`, resolved after the `set_knn` on `:1429` since `object` is already reassigned |
| `bbknn_sc` | `methods_sc_batch_corr.R:690, 717` | kNN gets `"pca"` (`embd_to_use` is `assertChoice(..., "pca")` at `:613`); sNN gets `"knn"` |
| `run_wnn_sc` | `methods_sc_multimodal.R:288` | kNN gets both source embedding names, modality-qualified; sNN gets `"wnn:knn"`. It builds `other_data[["wnn"]]` by hand, so stamp the two payloads there directly |
| `map_symphony_query` | `methods_sc_annotations_reference.R:557-559` | `character()`. The producer is a `SymphonyReference` which owns no cache, and it writes onto the **query**. `.stamp_artefact()` hashes the query's cell set, which is exactly right - this is the case that proves stamping belongs at the forwarder against the receiving object |

### Existing ad-hoc guards

- `run_paga_sc` (`methods_sc_trajectory.R:308-319`): keep. It compares a
  user-supplied obs column against the kNN, a different axis. Trim the advice
  from its message now that `assert_sc_state()` owns that.
- `assign_sc_type` (`methods_sc_annotations_reference.R:130-135`): keep, drop
  the "regenerate the neighbours" line for the same reason.
- `load_existing` ADT rows check (`methods_sc_io.R:1820`): unchanged, it guards
  `adt_counts`, which is not a cache artefact.

### `set_cells_to_keep`, `print`, options

`set_cells_to_keep` (`classes_single_cell_.R:1886`), after the DuckDB write:
warn once, listing `.sc_state_report(x)[stale == TRUE]` as
`artefact (modality)`. `SingleCellsMultiModal` does not override this method, so
it sweeps all four modality slots for free. The `load_existing` fallback
branches (`methods_sc_io.R:1706, 1846`) call the S3 `ScMap` method directly, not
this one, so loading never warns spuriously.

`print` for all four classes gets one `Stale artefacts:` line from the same
report, wrapped in `tryCatch()` - a print method must never fail.

### New: `reset_cells_to_keep(object, force = FALSE)`

The way back to a pristine object: every cell in the binary restored, cache
wiped. Generic in `R/base_generics_sc.R`, methods on `SingleCells` and
`SingleCellsMultiModal`. **No method** for `SingleCellsSubset` (its cell set
*is* the subset, resetting it would make it not a subset) or `MetaCells` (no
`cells_to_keep` concept); both get an explicit informative error rather than a
missing-method one.

Two facts checked rather than assumed, both of which make this cheap and safe:

- `SingleCellDuckDB$set_cells_to_keep` (`R/classes_single_cell_db.R:313`) runs
  `UPDATE obs SET to_keep = cell_idx IN (...)`. It is a boolean flag on every
  row, **not** a delete, and `get_cells_to_keep` reads `WHERE to_keep`. So a
  reset is fully recoverable and loses no obs data.
- `SingleCellCountData$get_shape()` returns `c(n_cells, n_genes)` from the
  binary header, and is what every ingest path already uses to set `dims`
  (`methods_sc_io.R:372` and eleven siblings). Use `get_shape()[1]` as the
  authority for "all cells found on the binary file" rather than the cached
  `dims` prop.

Flow:

1. `n_cells <- get_sc_rust_ptr(object)$get_shape()[1]`.
2. Collect what would be destroyed from `.sc_state_report(object)[present == TRUE]`,
   plus the HVG selection.
3. **Nothing cached and nothing filtered** - reset silently, return. Do not
   prompt for a no-op.
4. **`force = FALSE`**:
   - `!interactive()` - `stop()` telling the caller to pass `force = TRUE`.
     `readline()` returns `""` immediately under `Rscript`, `knitr` and
     `tinytest`, so prompting there would silently pick a branch. This is the
     one detail that makes the function safe to call from a script.
   - otherwise `readline()` with a prompt naming every artefact by modality and
     the HVG selection, defaulting to no. `readline()` is base, so no new
     dependency; `utils::menu()` would mean adding `utils` to `Imports`.
   - Declined - return the object **completely unchanged** and say so. No
     partial reset.
5. Reset:
   - `sc_map$cells_to_keep_idx <- NULL` rather than materialising
     `seq_len(n_cells) - 1`. `get_cells_to_keep.SingleCells:1247` already falls
     back to all cells on an empty slot, so this restores the genuine pristine
     state and costs nothing at 1M cells.
   - `sc_map$hvg_gene_indices <- NULL`. See the note below.
   - New R6 method `SingleCellDuckDB$reset_cells_to_keep()`, a single
     `UPDATE obs SET to_keep = TRUE`. Going through the existing
     `set_cells_to_keep()` would write a 1M-row temp table just to flip every
     flag back on.
   - `S7::prop(object, "sc_cache") <- new_sc_cache()`. The multimodal method
     also resets `adt_cache`, `atac_cache` and drops `other_data[["wnn"]]`.

**Judgement call, flag it if you disagree.** You said reset the `ScCache`, and
HVG lives in `ScMap`, not the cache. I am clearing it anyway: `find_hvg_sc`
reads `get_cells_to_keep()` at `methods_sc_processing.R:1158`, so an HVG
selection is every bit as cell-set dependent as the PCA, and leaving it behind
while calling the result pristine recreates exactly the inconsistency this whole
change exists to kill. It is named explicitly in the confirmation prompt, and it
is a one-line change if you want it kept.

Option `bixverse.cache_check`, three values: `"error"`, `"warn"`, `"none"`.
Default is unset, meaning the two-tier behaviour above. `"error"` promotes the
getter warnings to errors for anyone who wants the strict mode, `"none"`
silences everything including `assert_sc_state()`. No `.onLoad` write;
`getOption("bixverse.cache_check", NULL)` inline matches the repo, which has no
package options today. Document it in an `@section Options:` block on
`R/bixverse-package.R`.

### Backward compatibility

Nothing to do. `new_sc_cache()` is unchanged, `save_sc_exp_to_disk()` and
`load_existing()` are unchanged, and an artefact with no attribute reports
`stamped = FALSE, stale = FALSE`. An older bixverse reading a newer
`memory.qs2` gets payloads carrying an attribute it ignores.

No S7 property validator on `sc_cache`. A validator fires on every
`S7::prop<-`, including the four writes in the multimodal `load_existing` path
at `methods_sc_io.R:1803-1806`, and would hard-fail legacy loads.

## Companion change: `bixverse.gpu`

`~/repos/shared/bixverse.gpu` is a pure consumer. It never constructs an
`ScCache`, defines no method for any of these generics, never indexes the cache
positionally, and never serialises a `SingleCells` object. Since `ScCache` keeps
its exact structure under this design, the blast radius there is smaller still.
Two things to do, in a separate commit after this lands:

1. Add `from` to its six write-back sites so GPU artefacts join the chain
   instead of being recorded as roots: `R/sc_gpu.R:393` (`set_knn`, `from =
   embd_to_use`), `:422` (`set_snn_graph`, `from = "knn"`), `:529` / `:558`
   (same for the CAGRA path), `:780` (`set_embedding("harmony_gpu")`,
   `from = "pca"`), `:966` and `:1149` (UMAP / t-SNE, `from = embd_to_use`
   resolved against `cache_modality`, same rna-read/wnn-write asymmetry as the
   CPU path). PCA at `:674-676` stays a root. Without this, a GPU kNN is never
   marked stale when the PCA is re-run.
2. Consider `assert_sc_state()` at the entry of `find_neighbours_gpu_sc:339`,
   `find_neighbours_cagra_sc:484`, `.generate_seacells_gpu:191` and
   `.fast_cluster_gpu:199`, mirroring the CPU tier. Then bump the `bixverse`
   floor in its `DESCRIPTION` from `>= 0.4.0`; it carries `Remotes:
   GregorLueg/bixverse`, so its CI tracks bixverse HEAD.

Its tests should survive unchanged. The two that rely on forgiving getters,
`test_sc_gpu.R:126-133` (missing HVG) and `:255-262` (missing PCA), are
*absent*-artefact cases, and an absent artefact has no stamp. The sequences that
re-run `find_hvg_sc` or write an unrelated embedding mid-test (`:191-208`,
`:643-712`) are safe: nothing here invalidates on an HVG write, and writing
`"umap"` only mints a new id for `"umap"`. The two `bixverse:::` internals it
reaches into, `.get_manifoldsr_knn` and `.get_manifoldsr_knn_from_wnn`
(`R/sc_gpu.R:937, 939, 1121, 1123`), go through `get_knn_obj`, so the warning
tier lands in its UMAP/t-SNE paths automatically.

## Known limitations, stated up front

- An artefact without a stamp passes silently, so a mixed object stays partly
  unverified: re-run PCA on a legacy object and the new PCA is stamped while the
  pre-existing kNN stays `unknown`. `get_sc_cache_status()` reports it as
  `unknown` rather than pretending it is verified.
- The stamp records provenance, not parameters. A kNN built on the first 5 PCs
  and one built on all 15 of the same PCA carry the same `from`, so they are not
  distinguished. Bookkeeping gap, not a correctness hazard: no index mismatch
  can arise from it.
- Anything that pulls a payload out and rebuilds it loses the attribute. Inside
  the package that only happens on local copies after the read, but it means the
  stamp is not a durable provenance record for user code, only an internal
  consistency check. Stripping on the way out makes that explicit.
- Cluster labels live in the DuckDB obs table, not the cache, so they are not
  stamped. `run_paga_sc`'s hard error stays the guard for clusters versus kNN.

## Documentation

`devtools::document()`. Four new exports (`check_sc_state`, `assert_sc_state`,
`get_sc_cache_status`, `reset_cells_to_keep`), all of which must be added to
`_pkgdown.yml` in the *Single cell class and getters* block around line 336 -
pkgdown errors on an undocumented export. `reset_cells_to_keep` belongs next to
`set_cells_to_keep` there. `air format .` before committing. `R/zzz.R` needs no
change: stamps are plain lists with no class attribute, so no `print` method to
register. No `DESCRIPTION` change, `rlang` is already imported.

## Verification

1. `R CMD INSTALL --preclean .`, then `air format .`.
2. New `inst/tinytest/test_sc_cache_state.R`, following the setup in
   `inst/tinytest/test_sc_processing.R:1-70`:
   - Happy path `calculate_pca_sc -> umap_sc -> find_neighbours_sc`. Assert
     `get_sc_cache_status()` is a `data.table`, all `!stale`, and the chain links
     up: `umap$from == pca$id`, `knn$from == pca$id`, `snn$from == knn$id`. Ids
     unique and 16 chars.
   - **Attribute hygiene.** Assert `is.null(attributes(get_pca_factors(x))$bixverse_stamp)`
     and the same for `get_embedding`, `get_knn_obj`, `get_snn_graph`, so the
     stamp never leaks to user code.
   - **The reported scenario.** `set_cells_to_keep()` to a subset ->
     `expect_warning`, message names `pca`, `umap`, `knn`, `snn`. Assert
     `get_pca_factors()` warns and `find_neighbours_sc()` errors. Re-run
     `calculate_pca_sc` only -> PCA clean, `knn` and `snn` still stale,
     `generate_bt_meta_cells_sc()` still errors via the broken parent. Re-run
     `find_neighbours_sc` -> all fresh, `snn$from` points at the *new* `knn$id`.
   - Recursion: re-run PCA only, assert the sNN is stale *through* the kNN.
   - Options: `"error"` promotes the getter warning, `"none"` gives
     `expect_silent` even on the assert tier. Restore with `on.exit()`.
   - Legacy: strip the attributes off a populated cache by hand, assert all six
     getters still return and `get_sc_cache_status()` reports
     `stamped = FALSE, stale = FALSE`.
   - Round trip: `save_sc_exp_to_disk()` then `load_existing()`, assert ids
     survive byte-identically through both the `qs2` and `rds` paths.
   - `remove_knn()` leaves the sNN stale with no orphan stamp anywhere.
   - `expect_stdout(print(obj), "Stale artefacts")` for all four classes.
   - **`reset_cells_to_keep()`.** Tests run non-interactively, which is exactly
     the branch that must be nailed down: `expect_error(reset_cells_to_keep(obj))`
     on a populated object with `force = FALSE`, and `force = TRUE` succeeding.
     After a `force = TRUE` reset assert `length(get_cells_to_keep(obj)) ==
     get_sc_rust_ptr(obj)$get_shape()[1]`, that the DuckDB agrees
     (`get_sc_duckdb(obj)$get_cells_to_keep()` has the same length), that
     `get_sc_cache_status()` reports every artefact `present = FALSE`, and that
     `get_hvg()` warns and is empty. Then re-run the full
     `find_hvg_sc -> calculate_pca_sc -> find_neighbours_sc` chain on the reset
     object and assert it comes back all-fresh, proving the reset leaves a
     genuinely usable object rather than a half-wired one. Also assert the
     no-op path (nothing cached, nothing filtered) returns silently under
     `force = FALSE` without erroring. For multimodal, assert `adt_cache`,
     `atac_cache` and `other_data$wnn` are all cleared too.
3. Extend existing files rather than duplicating fixtures: `test_sc_subset.R`
   (fresh subset reports `present = FALSE`; running PCA on the subset leaves the
   parent's stamps untouched), `test_sc_aggregations.R` (MetaCells stays fresh
   after adding an obs column via `[[<-`), `test_sc_multimodal_adt.R`
   (`rna:pca` and `adt:pca` get different ids; narrowing cells warns for both
   modalities; `wnn:snn$from == wnn:knn$id` and `wnn:knn$from` has two entries).
4. Full `tinytest::test_package("bixverse")`. Every existing
   `set_cells_to_keep()` call site (`test_sc_processing.R:247`,
   `test_sc_seurat.R:184`, `test_sc_neighbours.R:75`, `test_sc_direct_load.R:47`,
   `test_sc_aggregations.R:547`, `test_sc_multimodal_adt.R:137, 151`) applies the
   filter before any PCA or kNN exists, so no new warning should fire.
   Re-confirm after implementation.
5. Install this branch, then run `tinytest::test_package("bixverse.gpu")` against
   it *before* touching the gpu package, to confirm the compatibility claim.
6. Render `vignettes/thinking_single_cell.qmd` and
   `vignettes/pbmc_single_cell.qmd`, the two that call `set_cells_to_keep()`
   mid-pipeline.
