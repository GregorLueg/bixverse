# Single cell: downstream analysis

Everything here assumes the core chain has run. Read `single-cell.md` first.
Most of these need a kNN, so `find_hvg_sc()` then `calculate_pca_sc()` then
`find_neighbours_sc()` is the entry price.

## Gene set scoring

| Method | Function | Params |
|---|---|---|
| module scores (Seurat-style) | `module_scores_sc()` | plain args |
| AUCell | `aucell_sc()` | `params_sc_aucell()` |
| VISION | `vision_sc()`, `vision_w_autocor_sc()` | `params_sc_vision()` |
| top expressed gene fraction | `top_genes_perc_sc()` | plain args |

`vision_w_autocor_sc()` adds spatial autocorrelation of the signature over the
kNN graph, which tells you whether a score is structured or noise.

Scores land as a data.table; push them onto the object with
`add_sc_new_obs(obj, obs_data = ...)`.

## Hotspot

Gene autocorrelation over the graph, then gene-gene local correlations, then
modules.

```r
obj <- hotspot_autocor_sc(obj, hotspot_params = params_sc_hotspot())
obj <- hotspot_gene_cor_sc(obj)
res <- generate_hotspot_membership(obj)
get_hotspot_membership(res)
```

## SCENIC and GRN inference

The full chain, with the CisTarget prerequisites:

```r
tf_list <- identify_tf_to_genes(...)
obj     <- scenic_gene_filter_sc(obj, ...)
grn     <- scenic_grn_sc(obj, scenic_params = params_scenic())
grn     <- tf_to_genes_correlations(grn, ...)

rankings <- read_motif_ranking(path)          # download_cistarget_hg38()
grn      <- tf_to_genes_motif_enrichment(grn, rankings,
                                         cistarget_params = params_cistarget())

regulons <- build_regulons(grn)
binary   <- binarise_regulon_activity(regulons,
                                      params = params_scenic_binarise())
```

Regressor choice lives in `params_scenic()`, with
`params_scenic_random_forest_defaults()`,
`params_scenic_gradient_boosting_defaults()` and
`params_scenic_extra_trees_defaults()` as the nested bundles.

Motif rankings are a large download. `download_cistarget_hg38()` fetches them.

## Differential abundance (miloR)

```r
milo <- get_miloR_abundances_sc(obj, milo_params = params_sc_miloR())
milo <- add_nhoods_info(milo, cell_info)
milo <- test_nhoods(milo, ...)

get_differential_abundance_res(milo)
get_index_cells(milo)
get_model_fit(milo)
```

The parameter bundle is `params_sc_miloR()`, with the capital R.
`generate_sc_knn()` builds a standalone kNN if you need one outside the object's
cache.

## Pseudobulk

```r
groups <- split(obs$cell_id, obs$sample_id)   # a NAMED list of cell ids
mat <- get_pseudobulked_sc(obj, cell_list = groups,
                           return_format = "dense", assay = "raw")
```

`cell_list` is a named list of cell identifiers, not a column name, so you build
the grouping yourself. `assay = "raw"` sums counts, `assay = "norm"` means the
normalised values. Returns aggregated cells x genes, ready for the bulk workflow
in `bulk.md`. This is the bridge between the two halves of the package.

## Trajectory inference

Order is fixed and each step asserts on the previous one:

```r
obj    <- run_magic_sc(obj, magic_params = params_sc_magic())
obj    <- run_paga_sc(obj, cluster_col = "leiden_clusters")
obj    <- run_palantir_sc(obj, palantir_params = params_sc_palantir())
trends <- run_gene_trends_sc(obj, gene_trends_params = params_sc_gene_trends())
```

`run_gene_trends_sc()` needs the MAGIC artefact, not just the kNN. Terminal
state selection is tuned with `params_sc_branch_selection()`.

`meld_sc()` (with `params_meld()`) scores per-cell condition density over the
graph, which is a different question from PAGA connectivity.

## NMF

```r
res <- nmf_sc(obj, nmf_params = params_nmf_hals())
res <- stabilised_nmf_sc(obj, ...)     # multi-run, picks a stable solution
get_w(res); get_h(res); get_best_run(res)
```

`get_hvg_data_sc()` gives you the HVG-restricted matrix if you want to factorise
that instead of the full thing.

## Reference mapping and annotation

**Symphony**, for projecting a query onto a built reference:

```r
ref   <- build_symphony_ref(reference_obj, ...)      # SymphonyReference
query <- map_symphony_query(query_obj, ref, symphony_params = params_symphony_map())
query <- transfer_labels_symphony(query, ref, ...)
query <- add_symphony_labels(query, ...)
```

**scType**, marker-based, no reference object needed:

```r
marker_list <- prepare_cell_markers(obj, cell_markers_dt)
scores      <- calc_sc_type_scores(obj, cell_marker_list = marker_list)
anno        <- score_clusters(scores, obj[[]][["leiden_clusters"]])
obj[["cell_type"]] <- anno$cell_type[match(clusters, anno$cluster_id)]
```

`cell_markers_dt` wants `cell_type` and `gene_id` columns, with `gene_id` in the
same identifier space as `var`. `assign_sc_type()` with `params_sctype_cells()`
does the per-cell variant.

## Ligand-receptor (NicheNet)

```r
influence <- generate_ligand_target_influence(networks,
                                              ligand_params = params_ligand_target())
lt_matrix <- get_influence(influence)
obj_info  <- compute_expression_info_sc(obj, grouping_variable = "cell_type")
activity  <- ligand_activity_scores(lt_matrix, geneset_of_interest, background)
top       <- prioritise_interactions(activity, obj_info, ...)
```

Networks come from outside bixverse. `genewalkR` is one source.

## Batch metrics

`calculate_kbet_sc()`, `calculate_batch_asw_sc()`, `calculate_batch_lisi_sc()`.
Run them on the uncorrected embedding first so you have a baseline, otherwise
the numbers mean nothing.

## Plot data extractors

bixverse computes, `bixverse.plots` draws. If you want the numbers rather than a
ggplot, or you don't have the plotting package installed:

`extract_embedding_data()`, `extract_dot_plot_data()`,
`extract_feature_plot_data()`, `extract_gene_violin_data()`,
`extract_gene_expression()`, `extract_paga_plot_data()`, `extract_feature_pair()`.

## Meta cells

Users never call `MetaCells()` directly. Three generators, all needing a kNN on
the parent:

```r
knn_object <- get_knn_obj(obj)              # needed for diffusion coordinates

mc <- generate_bt_meta_cells_sc(obj,        # hdWGCNA-style, bootstrapped
        sc_meta_cell_params = params_sc_bt_metacells(
          target_no_metacells = 250L, knn = list(k = 25L)),
        regenerate_knn = TRUE)

mc <- generate_seacells_sc(obj, params_sc_seacells())
mc <- generate_supercells_sc(obj, params_sc_supercell())
```

The important structural difference: unlike `SingleCells`, a `MetaCells` object
**holds its counts in memory**, with plain data.tables for obs and var and no
DuckDB. It shares `ScCache` and dispatches through the `ScOrMc` S7 union, so the
familiar generics work on it directly: `find_hvg_sc()`, `calculate_pca_sc()`,
`find_neighbours_sc()`, `find_clusters_sc()`, `umap_sc()`, `aucell_sc()`,
`vision_sc()`, `vision_w_autocor_sc()`.
`find_hvg_sc()` ignores `streaming` there and `calculate_pca_sc()` ignores
`sparse_svd`, since neither applies to an in-memory object.

Metacell-specific:

```r
mc  <- calc_diffusion_coordinates(mc, knn_data = knn_object)
mc  <- calc_meta_cell_purity(mc, grouping = "cell_type")
get_meta_cell_purity(mc)
calc_manifold_metrics(mc)                   # compactness, separation
merged <- merge_meta_cells(list(mc1, mc2))
mc_counts_to_list(mc); mc_get_clr_offsets(mc)
```

Merged objects set `is_merged = TRUE`, after which anything resolving back to
source cell indices bails out early. Pick which metacell method to use by
running the purity and manifold metrics on all three, that's what
`vignettes/meta_cells.qmd` does.

Metacells are also the natural input to the bulk co-expression module methods in
`bulk.md`, which is usually the point of making them.

## Multi-modal (CITE-seq)

`SingleCellsMultiModal` inherits from `SingleCells` and adds an ADT arm. The
order is forced:

```r
mm <- SingleCellsMultiModal(dir_data = dir)

mm <- load_tenx_h5(mm, ...)                    # RNA first
# ... QC, set_cells_to_keep() ...

adt <- read_tenx_h5_adt(path)                  # or read_multi_tenx_h5_adt()
adt <- detect_adt_isotypes(adt)
adt <- remove_adt_isotypes(adt)

mm <- add_adt_counts_sc(mm, adt, method = "clr")
# DSB instead, extra args pass through ... to new_adt_counts_dsb():
mm <- add_adt_counts_sc(mm, adt, method = "dsb",
                        dsb_params = params_sc_dsb(),
                        empty_drops = empty_drop_matrix,
                        isotype_names = isotypes)

# RNA arm
mm <- find_hvg_sc(mm); mm <- calculate_pca_sc(mm)
# ADT arm, note the separate generic
mm <- calculate_pca_adt_sc(mm, no_pcs = 10L)
mm <- find_neighbours_sc(mm, modality = "adt")

mm <- generate_wnn_graph_sc(mm, modality_1 = "rna", modality_2 = "adt",
                            wnn_params = params_sc_wnn())
mm <- find_clusters_sc(mm, modality = "wnn")
mm <- umap_sc(mm, modality = "wnn")
```

Constraints:

- ADT goes on **after** RNA ingest and after `cells_to_keep` is settled.
  `add_adt_counts_sc()` subsets to the kept cells and asserts they all exist in
  the ADT rownames.
- ADT normalisation (CLR or DSB) is baked into `add_adt_counts_sc()`, there is
  no separate step.
- Three caches, selected by the `modality` argument: `"rna"`, `"adt"`, `"wnn"`.
  Passing `modality != "rna"` to a generic on a plain `SingleCells` errors.
- `generate_wnn_graph_sc()` needs a PCA on **both** modalities.

Getters: `get_adt_names()`, `get_adt_feature_info()`, `get_adt_sample_info()`.
