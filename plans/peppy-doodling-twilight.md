# A usage skill for bixverse

## Context

bixverse has grown past the point where an agent can pick it up by reading the
code. The numbers:

- 654 `export()` entries in `NAMESPACE`, of which 179 are `rs_*` Rust wrappers
  users should not touch. Roughly 475 functions are genuinely user-facing.
- 62k lines of R across 70 files, ~22k of which is the single cell suite.
- 13 S7 classes, 5 R6 classes, 192 S7 generics, 366 method registrations.
- 61 `params_*()` constructors.
- 938 `.Rd` files, 23 vignettes totalling 11.7k lines.

An agent dropped into a user's analysis directory has none of this. It will
guess at signatures, invent arguments that do not exist, call things in the
wrong order, or reach for Seurat because that is what its training data is
thick with. The README already admits the problem: "there are methods hidden
here and there that lack any vignettes for now".

So yes, a skill makes sense, and this is close to the ideal case for one. The
package has a small number of strong, learnable conventions (`params_*()`
bundles, S7 chains that return the object, `get_results()`) wrapped around a
very large surface. That ratio is exactly what a skill is good at: teach the
conventions once, then route to a reference for the specifics.

Important: this is **not** `CLAUDE.md`. That file tells an agent how to
*develop* bixverse (cargo patches, `DEV_BUILD`, air, the extendr split). This
skill tells an agent how to *use* bixverse on someone's data, in some other
directory, with no repo checked out.

### Decisions taken

- **Lives in the package.** Authored at `inst/skills/bixverse/`, version
  controlled here, shipped on install. An exported `install_agent_skill()`
  copies it to `~/.claude/skills/bixverse/`. This is the only option that
  versions with the code and works outside the repo. `inst/extdata/*.parquet`
  already sets the precedent for shipping non-R assets, and `R/helpers.R:32`
  already has the `system.file()` lookup pattern to copy.
- **bixverse only.** The sister packages get a routing table (what they are,
  when to reach for them), not documented APIs. Five repos in sync is a
  maintenance trap.
- **Hybrid content.** Workflows, ordering constraints and gotchas are
  hand-written, because no generator can infer them. The flat function index is
  generated from `_pkgdown.yml` plus `.Rd` titles, because 654 exports rot fast.

## Structure

```
inst/skills/bixverse/
  SKILL.md                      # ~150 lines, always loaded on trigger
  references/
    install.md
    conventions.md
    single-cell.md
    single-cell-analysis.md
    enrichment.md
    bulk.md
    graphs-ontology.md
    api-index.md                # GENERATED, do not hand-edit
```

`SKILL.md` is the only thing that enters context on trigger, so it stays small.
Everything else is read on demand.

### SKILL.md

Frontmatter `name: bixverse`, and a `description` that fires on the package
name, on `SingleCells`/`MetaCells`/`BulkDge`/`BulkCoExp`, and on task phrasing
("single cell QC", "gene set enrichment in R", "co-expression modules").

Body, in order:

1. One paragraph on what the package is and the on-disk single cell model, so
   an agent does not assume Seurat-style in-memory objects.
2. Install one-liner plus the Windows caveat, pointing at `references/install.md`.
3. **The three invariants.** Tunables come in `params_*()` bundles, not long
   argument lists. Analysis classes are S7, every step returns the object so it
   pipes, results come out of `get_results()` / `get_outputs()`. Every exported
   function validates with checkmate, so a wrong argument errors loudly rather
   than silently doing the wrong thing.
4. **A router table**: user asks about X, read `references/Y.md`. This is the
   main job of the file.
5. A ten-line smoke test using `generate_single_cell_test_data()` so an agent
   can confirm the install works without downloading anything.
6. **Global traps**, short list: `rs_*` functions are internals with no
   validation, always use the R wrapper; the snake_case constructors in
   `R/deprecated.R` (`single_cell_exp`, `bulk_coexp`, `network_diffusions`, ...)
   are deprecated since 0.3.0 in favour of the UpperCamelCase S7 ones and older
   material still shows them; single cell normalisation happens at load, not as
   a step; a missing prerequisite warns and returns the object unchanged rather
   than erroring, so check the object after each step; `ScMap` indices are
   0-based for Rust; plotting mostly lives in `bixverse.plots`.
7. A pointer that the package is not Seurat or Scanpy and an agent should not
   translate idioms across. Where a Seurat equivalent exists the name is usually
   different and the semantics are not identical.

### references/install.md

The three README steps verbatim (Rust via rustup, `install.packages("rextendr")`,
`devtools::install_github`). Windows is unsupported, the blocker is h5
cross-compilation, do not attempt workarounds. Which `Suggests` packages are
needed for which method (Seurat, fgsea, GSVA, msigdbr, scDblFinder). The
ecosystem table: `bixverse.plots` for plotting, `bixverse.gpu` for GPU Harmony
/ PCA / kNN, `manifoldsR` for UMAP/tSNE/diffusion maps, `genewalkR` for
gene-gene networks.

### references/conventions.md

The cross-cutting grammar, which is what actually makes an agent competent here:

- The `params_*()` system: naming (`params_<method>`, `params_sc_<method>`,
  `params_*_defaults()` for nested zero-arg bundles), that they return validated
  named lists, and that they are passed as one argument. Include the nesting
  quirk, because it is not guessable: `params_sc_neighbours(knn = list(k = 20L))`
  takes the `knn` list, `modifyList`s it over `params_knn_defaults()` and then
  **flattens** it, so downstream code reads `params$k`, not `params$knn$k`.
- Which steps take a bundle and which take plain arguments. It is not uniform:
  `find_clusters_sc()`, `find_markers_sc()`, `umap_sc()`, `run_cell_qc()`,
  `set_cells_to_keep()` and `get_pseudobulked_sc()` all take plain arguments,
  while nearly everything else takes a bundle. A table is the right format.
- Bundles are validated by matching `assert*` functions in
  `R/checkmate_extensions.R`, so a malformed bundle fails at the call site with
  a useful message. Worth saying, it tells an agent that guessing is safe-ish.
- The S7 chain shape: `ClassConstructor()` then `preprocess_*()` then verb
  methods then `get_results()`.
- Verb vocabulary: `get_` never mutates, `calculate_`/`calc_` computes and
  stores back, `find_` discovers, `run_` runs a named algorithm, `add_` appends.
- `[[` and `[[<-` on `SingleCells` are sugar for obs columns.
- Every class has a `print` method, so printing the object is the cheapest way
  to see state. Tell the agent to do this rather than guess.
- `BixverseBaseClass` and the `get_params()` / `get_results()` contract.

### references/single-cell.md

The highest-value file in the set. The core workflow in call order, then the
traps, which is where an agent actually goes wrong.

```
SingleCells(dir_data)              # one arg, the dir must already exist
  -> load_mtx() | load_h5ad() | load_tenx_h5() | load_seurat() | load_r_data()
     | stream_h5ad()   # multi-sample: prescan_*() -> load_multi_*()
  -> gene_set_proportions_sc() / top_genes_perc_sc()
  -> scrublet_sc() | scdblfinder_sc() | doublet_detection_boost_sc()
  -> run_cell_qc() -> set_cells_to_keep()
  -> find_hvg_sc() | find_hvg_batch_aware_sc()
  -> calculate_pca_sc()
  -> [harmony_sc | harmony_v2_sc | fast_mnn_sc | bbknn_sc | seurat_cca_sc
      | seurat_rpca_sc]            # writes a NEW named embedding
  -> find_neighbours_sc(embd_to_use = "pca" | "harmony" | ...)
  -> find_clusters_sc() | fast_cluster_sc()
  -> umap_sc() | tsne_sc() | phate_sc()
  -> find_markers_sc() | find_all_markers_sc() | find_specific_markers_sc()
  -> save_sc_exp_to_disk()
```

The traps, in rough order of how much damage they do:

- **There is no `normalise_sc()`.** Log-CPM normalisation runs inside Rust
  during `load_*`, driven by `target_size` in `params_sc_min_quality()`. The
  `min_lib_size` / `min_unique_genes` / `min_cells` cutoffs are applied at load
  too and are **not reversible without re-loading**, because failing cells and
  genes are never written to the binary. An agent looking for a normalisation
  step will not find one and may invent one.
- **Missing prerequisites warn and return the object unchanged.**
  `calculate_pca_sc()` without an HVG selection is a silent no-op
  (`R/methods_sc_processing.R:1268`), same for `find_neighbours_sc()` without
  the named embedding (`:1404`). The pipeline appears to run and produces
  nothing. Table of every step and its hard prerequisite goes here.
- **Mixed mutation semantics.** `setnames_sc()` and `drop_cols_sc()` mutate the
  DuckDB in place and return invisibly, so they are called without
  reassignment. Everything else mutates the in-memory S7 object and **must** be
  reassigned. Both styles appear in the same workflow.
- **Filter first, then HVG and PCA.** `set_cells_to_keep()` after a PCA does not
  delete the PCA, it stamps it stale, and the next `find_neighbours_sc()`
  errors. `R/functions_sc_cache_state.R` hashes the cell set and walks parent
  artefact ids, so staleness is transitive. Inspect with
  `get_sc_cache_status()`; the tier is controlled by
  `options(bixverse.cache_check = "error" | "warn" | "none")`. Tell the agent to
  read the warning and recompute from the stale step, never to silence it.
- **0-based versus 1-based.** `get_cells_to_keep()` and `get_hvg()` return
  0-based Rust indices. `set_cells_to_keep()` / `set_hvg()` accept characters or
  1-based numerics and subtract 1 (`R/classes_single_cell_.R:99-131`).
  `calculate_pca_sc(hvg = ...)` expects 1-based. This asymmetry is the single
  easiest thing for an agent to get wrong.
- **`obj[["col"]] <- x` matches by position, not name**, and requires
  `length(x) == length(get_cells_to_keep(obj))`. Wrong length errors; wrong
  order silently mislabels every cell.
- **`obj[["col"]]` is filtered, `get_sc_obs(obj)` is not** (`filtered = FALSE`
  by default). Mixing them when building QC vectors gives length mismatches.
- `reset_cells_to_keep()` wipes the cache and the HVG selection, and errors in
  non-interactive sessions without `force = TRUE`.
- The cache is memory-only. `load_existing()` restores counts and obs/var, but
  PCA/kNN/embeddings/HVG come back only if `save_sc_exp_to_disk()` ran first.
- Reserved embedding names: `pca`, `knn`, `snn`, `magic`. Do not use them as
  `slot_name`.
- `streaming` is an integer (`0L`/`1L`/`2L`) at ingest but a
  `NULL`/logical on `find_hvg_sc()`, `aucell_sc()` and friends. Same name,
  different meaning.
- Silent auto-switches: `calculate_pca_sc()` flips to `sparse_svd` above 500k
  cells, `find_all_markers_sc()` downsamples the rest-group to 100k unless
  `downsampling = FALSE`.
- **On-disk, not in memory.** `dir_data` must persist for the object to stay
  usable, so `tempdir()` is for experiments only.

Then the pipeline DSL:

```r
pipe <- sc_pipeline() %>>% step_hvg_sc() %>>% step_pca_sc() %>>%
  step_harmony_sc() %>>% step_neighbours_sc() %>>% step_clusters_sc()
apply_pipeline(pipe, sc_obj)                          # validates, then runs
apply_pipeline_per_group(pipe, sc_obj, group_col = "sample_id")
```

Be precise about what `validate_pipeline()` (`R/pipelines_sc.R:192`) does, since
overselling it would be worse than not mentioning it: `apply_pipeline()` calls
it up front and it checks **class dispatch** (can each step run on what the
previous step returns), naming the offending step. It does **not** check
step ordering, so an HVG-after-PCA pipeline passes validation and then hits the
silent no-op above. Nine steps exist: `step_hvg_sc`, `step_pca_sc`,
`step_neighbours_sc`, `step_clusters_sc`, `step_harmony_sc`,
`step_harmony_v2_sc`, `step_bbknn_sc`, `step_fast_mnn_sc`, `step_metacells_sc`.

Then `SingleCellsSubset`: `SingleCellsSubset(sc_object, grouping_column, group)`
shares the parent's count pointer, keeps indices in the **parent's** index
space, gets a fresh empty cache and **drops the parent's HVG** deliberately
(`R/classes_single_cell_subset.R:99`), so the chain restarts at `find_hvg_sc()`.
`merge_subset_obs()` pushes labels back up.

Finally the writers: `save_sc_exp_to_disk()`, `merge_sc_experiments()`,
`write_h5ad_sc()`, `write_cellranger_output()`.

A minimal working example closes the file, built on
`generate_single_cell_test_data()` so it needs no network. The shape is in
`inst/tinytest/helper_sc.R` (`sc_test_fixture()`, `sc_test_object()`,
`sc_test_prepped()`); the reference version must inline the real calls rather
than depend on those internal test helpers.

### references/single-cell-analysis.md

The downstream branches, each as a short chain: gene set scoring (`aucell_sc`,
`module_scores_sc`, `vision_sc`, `top_genes_perc_sc`), Hotspot, SCENIC/GRN
(including the CisTarget prerequisites), miloR differential abundance,
pseudobulk, trajectory (MAGIC -> PAGA -> Palantir -> gene trends), single cell
NMF, reference mapping (Symphony and scType), NicheNet, batch metrics (kBET,
ASW, LISI), and the `extract_*_data()` family that feeds `bixverse.plots`.

Most of these need a kNN, so the HVG -> PCA -> neighbours chain has to have run:
`run_paga_sc`, `run_palantir_sc`, `meld_sc`, `get_miloR_abundances_sc` and all
three metacell generators assert on it. `run_gene_trends_sc` additionally needs
`run_magic_sc`.

**Meta cells.** Users never call `MetaCells()` directly, they get one from
`generate_bt_meta_cells_sc()`, `generate_seacells_sc()` or
`generate_supercells_sc()`. The key difference to state plainly: unlike
`SingleCells`, a `MetaCells` object **holds its counts in memory**, with plain
data.tables for obs/var and no DuckDB. It shares `ScCache` and dispatches
through the `ScOrMc` union in `R/classes_union.R`, so `find_hvg_sc()`,
`calculate_pca_sc()`, `find_neighbours_sc()`, `find_clusters_sc()`, `umap_sc()`
and `aucell_sc()` all work on it. Metacell-only:
`calc_meta_cell_purity()`, `calc_diffusion_coordinates()`,
`calc_manifold_metrics()`, `merge_meta_cells()`. Merged objects set
`is_merged = TRUE`, after which anything resolving back to source cell indices
bails out.

**Multi-modal.** `SingleCellsMultiModal(dir_data)` inherits from `SingleCells`.
Order is forced: ingest RNA first, settle `cells_to_keep`, then
`add_adt_counts_sc()`, which asserts the ADT rownames cover the kept cells and
which bakes in the normalisation (CLR or DSB, no separate step). ADT gets its
own generic `calculate_pca_adt_sc()`, and there are three caches selected by the
`modality` argument (`rna` / `adt` / `wnn`). `generate_wnn_graph_sc()` needs a
PCA on both modalities. Passing `modality != "rna"` to a generic on a plain
`SingleCells` errors.

### references/enrichment.md

Hypergeometric (`gse_hypergeometric*`, `simplify_hypergeom_res`), the GO
DAG-aware elimination path (`GeneOntologyElim` and `gse_go_elim_method*`,
`fgsea_go_elim`), rank-based (`calc_fgsea`, `calc_gsea_traditional`,
`calc_mitch`), per-sample scoring (`calc_gsva`, `calc_ssgsea`, `calc_singscore*`),
and where the GO data comes from (`get_go_data_human()` /
`load_go_human_data()`, backed by `inst/extdata/*.parquet`).

### references/bulk.md

`BulkDge` end to end (QC, preprocess, normalise, batch correct, PCA, limma /
limma-voom / Hedges' g, getters). `BulkCoExp` with the five module families
(CoReMo, graph Leiden, ICA, NMF/HALS, DGRDL) as parallel chains, plus
differential correlation and contrastive PCA. Guidance on which family to pick,
since the vignette benchmarks them against planted modules.

### references/graphs-ontology.md

`NetworkDiffusions` (diffuse, tied diffusion, permutation, AUC, community
detection), `RbhGraph`, `SimilarityNetworkFusion`, and both ontology paths (flat
functions on a parent-child table, and the `OntologySim` class).

### references/api-index.md (generated)

Flat index, grouped by the `_pkgdown.yml` reference sections, one line per
function: `- fn_name(): <Rd title>`. Roughly 560 lines once the "Rust wrappers"
section is collapsed to a single note. This is what an agent greps when it needs
to know whether a thing exists.

Note that `_pkgdown.yml` indexes 92 `rs_*` functions while `NAMESPACE` exports
179 of them, so the generator has to handle exports with no pkgdown home rather
than assume full coverage. Pushing those into an "unindexed" bucket makes the
gap visible instead of silently dropping them.

## Source material

Nothing here needs inventing. Write the references from:

- `vignettes/thinking_single_cell.qmd` and `vignettes/design_single_cell.qmd`,
  which are already the conceptual primers, for the mental model in
  `single-cell.md`.
- `vignettes/pbmc_single_cell.qmd` for the canonical chain, plus the 13 other
  single cell vignettes for the branches.
- `vignettes/gse_methods.qmd`, `pathway_activity.qmd`,
  `bulk_coexpression_modules.qmd`, `genetic_diffusions.qmd`, `ontologies.qmd`,
  `cpca.qmd` for the non-single-cell references.
- `inst/tinytest/` (47 files) for the minimal working examples. These are
  executable and current, which vignette snippets are not always. Prefer them
  for anything that goes into the skill as a runnable example.
- `R/functions_sc_cache_state.R` and `R/pipelines_sc.R` for the two mechanisms
  worth explaining in detail rather than listing.

The references are a compression of ~11.7k lines of vignette into something an
agent reads in one go, not a rewrite. Keep them dense: chains, tables, traps.
No prose warm-up.

## Files

New:

- `inst/skills/bixverse/SKILL.md` and `inst/skills/bixverse/references/*.md`
- `R/functions_agent_skill.R`, one exported function:

  ```r
  install_agent_skill(dest = "~/.claude/skills", overwrite = FALSE)
  ```

  checkmate-asserts both args, resolves the source with
  `system.file("skills", "bixverse", package = "bixverse")`, errors with a clear
  message if that is `""` (same shape as `R/helpers.R:37`), refuses to clobber
  an existing target unless `overwrite = TRUE`, copies recursively, messages the
  destination, returns the path invisibly. Roxygen with `@examples` wrapped in
  `\dontrun{}`.

- `data-raw/generate_api_index.R`, which reads `_pkgdown.yml`
  (`yaml::read_yaml`, dev-only dependency), pulls titles from
  `tools::Rd_db("bixverse")` and writes
  `inst/skills/bixverse/references/api-index.md`. It warns on any export missing
  from `_pkgdown.yml`, which doubles as a pkgdown coverage check.

Modified:

- `_pkgdown.yml`: add `install_agent_skill` to the "Utils" section, otherwise
  `R CMD check` flags an undocumented topic in the reference index.
- `NAMESPACE` and `man/`, regenerated by `devtools::document()`, not
  hand-edited.

Not modified: `.Rbuildignore` stays as is. `inst/` ships by design, and the
skill is markdown, so the size cost is negligible.

## Verification

1. `Sys.setenv(DEV_BUILD = "true")`, `devtools::document()`, then
   `R CMD INSTALL --preclean .`.
2. `system.file("skills", "bixverse", package = "bixverse")` returns a non-empty
   path containing `SKILL.md` and `references/`.
3. `bixverse::install_agent_skill(dest = tempdir())` copies the tree; a second
   call errors; with `overwrite = TRUE` it succeeds.
4. `Rscript data-raw/generate_api_index.R` on a clean tree leaves
   `git diff --stat` empty, proving the checked-in index matches the source.
5. Run the `SKILL.md` smoke snippet verbatim in a fresh R session. It uses
   `generate_single_cell_test_data()`, so it needs no network.
6. Cross-check every function name in the hand-written references against
   `NAMESPACE`. A skill that names a function that does not exist is worse than
   no skill. Do this with a script, not by eye.
7. End-to-end: install the skill, start Claude Code in an empty scratch
   directory (not this repo, so `CLAUDE.md` cannot leak in), and ask for
   something like "load a 10x mtx directory with bixverse, run QC, cluster it".
   Confirm the skill triggers, the reference gets read, and the emitted code
   runs. Repeat for a non-single-cell task (GO enrichment on a gene list) to
   check the router table works.
8. `air format .` and the pre-commit hooks on the new R files.

## Maintenance

The generated index is the only part that must be regenerated. Add
`Rscript data-raw/generate_api_index.R` to the release checklist alongside the
`DESCRIPTION` and `Cargo.toml` version bump, so it cannot drift silently across
a release.
