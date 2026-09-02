# Single cell: the core workflow

## The mental model

A `SingleCells` object holds no counts in R memory. It's a handle over a
directory:

- `counts_cells.bin`, cell-major (CSR), for per-cell operations
- `counts_genes.bin`, gene-major (CSC), for per-gene operations
- `sc_duckdb.db`, the obs and var tables

Plus two lightweight side-cars living on the S7 object itself: `ScMap` (name to
index maps, cells-to-keep, HVG indices) and `ScCache` (PCA factors, embeddings,
kNN, sNN graph).

Consequences you have to work with:

- `dir_data` must survive as long as the object. `tempdir()` is for experiments.
- Counts are never fully materialised in R. Pull what you need with
  `get_sc_counts()` or the `extract_*` helpers.
- The cache is memory-only. It survives a reload only if you called
  `save_sc_exp_to_disk()` first.

## The chain

```r
obj <- SingleCells(dir_data = some_existing_dir)   # one argument, dir must exist

# 1. ingest, exactly one loader
obj <- load_mtx(obj, sc_mtx_io_param = get_cell_ranger_params(path))

# 2. QC metrics onto the obs table
obj <- gene_set_proportions_sc(obj, gs_of_interest)
obj <- top_genes_perc_sc(obj, top_n_vals = c(5L, 10L, 25L))

# 3. doublets, optional
obj <- scrublet_sc(obj, params_scrublet())

# 4. call outliers and filter
qc <- run_cell_qc(metrics, cells_to_keep = get_cells_to_keep(obj), directions,
                  threshold = 3)
obj[["outlier"]] <- qc$combined
obj <- set_cells_to_keep(obj, cell_ids_of_the_keepers)

# 5. feature selection and reduction
obj <- find_hvg_sc(obj, hvg_no = 2000L)
obj <- calculate_pca_sc(obj, no_pcs = 30L)

# 6. batch correction, optional, writes a NEW named embedding
obj <- harmony_sc(obj, batch_column = "sample_id")

# 7. graph and clusters
obj <- find_neighbours_sc(obj, embd_to_use = "harmony")
obj <- find_clusters_sc(obj, res = 1, name = "leiden_clusters")

# 8. embedding for plots
obj <- umap_sc(obj)

# 9. markers
markers <- find_all_markers_sc(obj, column_of_interest = "leiden_clusters")

# 10. persist the in-memory cache
save_sc_exp_to_disk(obj, type = "qs2")
```

## Loading

| Input | Function |
|---|---|
| Cell Ranger MTX directory | `load_mtx()`, with `get_cell_ranger_params(path)` |
| h5ad | `load_h5ad()`, or `load_h5ad_norm()` for pre-normalised, or `stream_h5ad()` for big files |
| 10x HDF5 | `load_tenx_h5()` |
| Seurat object | `load_seurat()` |
| `SingleCellExperiment` | `load_sce()`, with `assay_name` naming the raw counts |
| in-memory sparse matrix + obs + var | `load_r_data()` |
| a directory written earlier | `load_existing()` |

Multi-sample variants: `load_multi_h5ad()`, `load_multi_mtx()`,
`load_multi_tenx_h5()`, with `prescan_h5ad_files()`, `prescan_mtx_dirs()` and
`prescan_tenx_h5_files()` to work out the shape first.

All loaders take `sc_qc_param = params_sc_min_quality()`:

```r
params_sc_min_quality(
  min_unique_genes = 100L,   # per cell
  min_lib_size     = 250L,   # per cell
  min_cells        = 10L,    # per gene
  target_size      = 1e4     # normalisation target
)
```

**Normalisation happens here.** Log-CPM to `target_size` runs inside Rust during
the load. There is no `normalise_sc()` and you should not write one.

`load_sce()` takes `colData` as obs and `rowData` as var. It does not carry
`reducedDims` or `altExps` over, and the assay you point it at has to be raw
counts. A lot of objects in the wild only ship `logcounts`, so name the assay
rather than trusting the default.

**The load-time cutoffs are irreversible.** Cells and genes failing
`min_unique_genes` / `min_lib_size` / `min_cells` are never written to the
binary. Getting them back means re-loading. Set them permissively and do the
real QC afterwards with `run_cell_qc()`, which is reversible.

`streaming` at ingest is an integer: `0L` in-memory, `1L` cell-batched
(default), `2L` gene-capped. Note that the same argument name on `find_hvg_sc()`
and `aucell_sc()` is a logical and means something different.

## QC

```r
qc_df <- obj[[c("cell_id", "lib_size", "nnz", "MT")]]

metrics <- list(
  log10_lib_size = log10(qc_df$lib_size),
  log10_nnz      = log10(qc_df$nnz),
  MT             = qc_df$MT
)
directions <- c(log10_lib_size = "twosided", log10_nnz = "twosided",
                MT = "above")

qc <- run_cell_qc(metrics, cells_to_keep = get_cells_to_keep(obj),
                  directions = directions, threshold = 3)
```

`run_cell_qc()` returns a `CellQc` S3 object carrying the metrics, per-metric
outlier calls and `$combined`. `run_cell_qc_fixed()` is the hard-threshold
version. `violin_plot_sc(qc)` and `joint_plot_sc(qc)` plot it, both from
`bixverse.plots`.

Then filter:

```r
obj[["outlier"]] <- qc$combined
obj <- set_cells_to_keep(obj, qc_df[!qc$combined, cell_id])
```

`cells_to_keep` is a mask, not a deletion. Counts stay on disk and every
downstream method respects it. `reset_cells_to_keep()` undoes it, but it wipes
the cache and the HVG selection, and errors in a non-interactive session unless
`force = TRUE`.

## Traps

### Missing prerequisites warn, they do not error

`calculate_pca_sc()` with no HVG selection warns "Returning object as is" and
hands the object back untouched. Same for `find_neighbours_sc()` without the
named embedding. A chain that looks like it ran can have done nothing. Print the
object or call `get_sc_cache_status()` after each step you're unsure of.

| Step | Needs |
|---|---|
| anything | a `load_*` call |
| `calculate_pca_sc` | `find_hvg_sc`, or explicit `hvg =` |
| batch correction | `calculate_pca_sc` |
| `find_neighbours_sc` | an embedding matching `embd_to_use` |
| `find_clusters_sc` | `find_neighbours_sc` (the sNN) |
| `umap_sc(use_knn = TRUE)` | a kNN |
| metacell generators, `run_paga_sc`, `run_palantir_sc`, `meld_sc`, miloR | a kNN |
| `run_gene_trends_sc` | `run_magic_sc` |

### Filter first, then HVG and PCA

Every cache artefact carries a provenance stamp recording the cell-set hash and
its parent artefact ids. Call `set_cells_to_keep()` after computing a PCA and
the PCA isn't deleted, it's stamped stale, transitively taking the kNN and sNN
with it. The next `find_neighbours_sc()` then errors.

```r
get_sc_cache_status(obj)   # one row per artefact: stamped, stale, reason
```

Two tiers by default: getters warn, `assert_sc_state()` errors.
`options(bixverse.cache_check = "error" | "warn" | "none")` overrides. Don't set
`"none"` to make a warning go away, recompute from the step that went stale.

### 0-based versus 1-based

`ScMap` stores Rust indices.

| Call | Indexing |
|---|---|
| `get_cells_to_keep()` | returns **0-based** |
| `get_hvg()` | returns **0-based** |
| `set_cells_to_keep()` | takes cell ids, or **1-based** numerics |
| `set_hvg()` | takes gene ids, or **1-based** numerics |
| `calculate_pca_sc(hvg = ...)` | takes **1-based** |

Prefer character identifiers over indices wherever a function accepts them. It
sidesteps the whole problem.

### Mixed mutation semantics

`setnames_sc()` and `drop_cols_sc()` mutate the DuckDB in place and return
invisibly, so they're called *without* reassignment:

```r
setnames_sc(obj, table = "var", old = "column1", new = "gene_symbol")
```

Everything else mutates the in-memory S7 object and **must** be reassigned:

```r
obj <- find_hvg_sc(obj, hvg_no = 2000L)
```

Both appear in the same workflow. Reassigning a `setnames_sc()` is harmless;
forgetting to reassign a `find_hvg_sc()` silently loses the work.

### Obs column assignment is positional

```r
obj[["cell_type"]] <- vec
```

requires `length(vec) == length(get_cells_to_keep(obj))` and matches by
**position**, not by name. Wrong length errors. Wrong order silently mislabels
every cell, which is much worse. Build the vector from
`get_sc_obs(obj, filtered = TRUE)` so the order is guaranteed.

Related: `obj[[...]]` is filtered, `get_sc_obs(obj)` defaults to
`filtered = FALSE`. Mixing them gives length mismatches.

### Reserved embedding names

`pca`, `knn`, `snn`, `magic`. Don't use them as a `slot_name`.

### Silent auto-switches

`calculate_pca_sc()` flips to sparse SVD above 500k cells to keep memory down.
`find_all_markers_sc()` downsamples the rest-group to 100k cells unless
`downsampling = FALSE`. Both message, neither errors.

## Clustering and markers

```r
obj <- find_clusters_sc(obj, cluster_algorithm = "leiden", res = 1,
                        name = "leiden_clusters")
```

For large data, `fast_cluster_sc()` runs mini-batch k-means, then kNN on the
centroids, then Louvain across a resolution grid. Optionally across seeds for
stability:

```r
res <- fast_cluster_sc(obj, resolutions = c(5, 3, 2, 1.5, 1, 0.5),
                       return_kmeans = TRUE, no_seeds = 25L, grid_search = TRUE)
obj <- add_sc_new_obs(obj, obs_data = get_data(res))
```

Markers: `find_all_markers_sc()` (each cluster vs rest),
`find_markers_sc(cells_1, cells_2)` (an explicit contrast),
`find_specific_markers_sc()` (genes specific to one group).

## Batch correction

All of these run after `calculate_pca_sc()` and write a new named embedding,
which you then hand to `find_neighbours_sc(embd_to_use = ...)`.

| Function | Embedding | Params |
|---|---|---|
| `harmony_sc()` | `"harmony"` | `params_sc_harmony()` |
| `harmony_v2_sc()` | `"harmony"` | `params_sc_harmony_v2()` |
| `fast_mnn_sc()` | `"mnn"` | `params_sc_fastmnn()` |
| `seurat_cca_sc()` | | `params_sc_seurat_cca()` |
| `seurat_rpca_sc()` | | `params_sc_seurat_rpca()` |
| `bbknn_sc()` | modifies the graph directly | `params_sc_bbknn()` |

Quantify the damage before and after with `calculate_kbet_sc()`,
`calculate_batch_asw_sc()` and `calculate_batch_lisi_sc()`.

## Pipelines

A declarative alternative to the call chain, useful when the same recipe runs
over many groups:

```r
pipe <- sc_pipeline() %>>%
  step_hvg_sc(hvg_no = 2000L) %>>%
  step_pca_sc(no_pcs = 30L) %>>%
  step_harmony_sc(batch_column = "sample_id") %>>%
  step_neighbours_sc() %>>%
  step_clusters_sc(name = "leiden_clusters")

obj <- apply_pipeline(pipe, obj)
per_group <- apply_pipeline_per_group(pipe, obj, group_col = "sample_id")
```

Nine steps: `step_hvg_sc`, `step_pca_sc`, `step_neighbours_sc`,
`step_clusters_sc`, `step_harmony_sc`, `step_harmony_v2_sc`, `step_bbknn_sc`,
`step_fast_mnn_sc`, `step_metacells_sc`.

`apply_pipeline()` runs `validate_pipeline()` up front. Be clear on what that
buys you: it checks **class dispatch**, that each step has a method for whatever
the previous step returns, and names the offending step if not. It does **not**
check step ordering, so an HVG-after-PCA pipeline validates fine and then hits
the silent no-op described above. Get the order right yourself.

Pipelines are inert until applied. Inspect one with `print(pipe)` or
`pipe$steps`.

## Subsets

For sub-clustering a cell type without touching the parent:

```r
sub <- SingleCellsSubset(obj, grouping_column = "cell_type", group = "T cells")
```

It shares the parent's count pointer (no copy) and keeps indices in the
**parent's** index space, so no translation is needed. It gets a fresh empty
cache and **deliberately drops the parent's HVG**, because highly variable genes
within one cell type are not the same set. So the chain restarts at
`find_hvg_sc()`.

Push labels back up with `merge_subset_obs()`.

## Writing out

`save_sc_exp_to_disk()` persists the in-memory cache (`type = "qs2"` is faster
and needs the `qs2` package). `merge_sc_experiments()` combines objects.
`write_h5ad_sc()`, `write_h5ad_sc_dense()`, `write_tenx_h5_sc()` and
`write_cellranger_output()` export to the usual interchange formats.

## Minimal working example

No downloads, runs in seconds.

```r
library(bixverse)

dir_sc <- file.path(tempdir(), "bixverse_mwe")
dir.create(dir_sc, showWarnings = FALSE)

syn <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 1000L, n_genes = 100L),
  seed = 42L
)

obj <- load_r_data(
  object = SingleCells(dir_data = dir_sc),
  counts = syn$counts,   # dgRMatrix, cells x genes
  obs = syn$obs,
  var = syn$var,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 45L, min_lib_size = 300L, min_cells = 500L
  ),
  streaming = 0L
)

obj <- find_hvg_sc(obj, hvg_no = 30L)
obj <- calculate_pca_sc(obj, no_pcs = 10L)
obj <- find_neighbours_sc(
  obj,
  neighbours_params = params_sc_neighbours(knn = list(knn_method = "exhaustive"))
)
obj <- find_clusters_sc(obj, res = 1, name = "clusters")

markers <- find_all_markers_sc(obj, column_of_interest = "clusters")

unlink(dir_sc, recursive = TRUE, force = TRUE)
```

For real data without a download of your own: `download_pbmc3k()`,
`download_pbmc8k()`, `download_pbmc_batches()`, `download_cd34_data()`,
`download_marrow_cd34()`, `download_demuxlet_pbmc()`,
`download_pbmc_totalseq_data()`.
