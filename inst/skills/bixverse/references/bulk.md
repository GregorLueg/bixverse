# Bulk RNAseq

Two classes, both S7, both following the standard chain. `BulkDge` for
differential expression, `BulkCoExp` for co-expression modules.

Pseudobulked single cell data feeds straight into either, see
`get_pseudobulked_sc()` in `single-cell-analysis.md`.

## Differential expression: BulkDge

```r
obj <- BulkDge(...)                    # or bulk_dge_from_h5ad(path, ...)

# metadata surgery, all optional
obj <- add_new_metadata(obj, ...)
obj <- update_metadata_values(obj, ...)
obj <- fix_meta_data_column(obj, ...)
obj <- change_gene_identifier(obj, ...)
obj <- remove_samples(obj, ...)

obj <- qc_bulk_dge(obj)
obj <- preprocess_bulk_dge(obj)
obj <- normalise_bulk_dge(obj)
obj <- batch_correction_bulk_dge(obj, batch_column = ...)
obj <- calculate_pca_bulk_dge(obj)

obj <- calculate_dge_limma(obj, ...)   # limma
obj <- calculate_dge_hedges(obj, ...)  # Hedges' g effect sizes
obj <- calculate_all_dges(obj, ...)    # everything across all contrasts
```

`run_limma_voom()` is the standalone voom path and is what
`calculate_dge_limma()` calls internally for count data.

Getters: `get_dge_list()`, `get_dge_limma_voom()`, `get_dge_effect_sizes()`,
`get_model_fit()`, `get_tpm_counts()`, `get_fpkm_counts()`,
`get_dge_qc_plot()`.

Plots: `plot_pca_res()` on the object. Several other `plot_*` helpers exist in
the source but are not exported, so reach for `get_dge_qc_plot()` and the
getters above and plot the returned data yourself.

## Co-expression modules: BulkCoExp

One class, five algorithm families, all starting the same way:

```r
coexp <- BulkCoExp(raw_data = data_log, meta_data = meta_data)
coexp <- preprocess_bulk_coexp(coexp, hvg = 0.5)
```

`hvg = 0.5` keeps the top half of genes by MAD. Note that **the object mutates
as each method runs**, so if you want to compare families on the same data, make
a copy per branch rather than chaining them off one object.

### Which family?

| Family | Good for | Cost |
|---|---|---|
| CoReMo | correlation-driven modules with a stability check | needs an epsilon choice, see below |
| Graph / Leiden | modules as graph communities, sub-clustering built in | resolution sweep |
| ICA | overlapping modules, latent factors | slowest, needs a component sweep |
| NMF (HALS) | parts-based, non-negative, interpretable loadings | needs non-negative input |
| DGRDL | sparse dictionary learning | grid search over two parameters |

`vignettes/bulk_coexpression_modules.qmd` benchmarks all five against planted
modules on synthetic data, which is the honest way to pick.

### CoReMo

```r
coexp <- cor_module_processing(coexp, cor_method = "spearman")
coexp <- cor_module_check_epsilon(coexp, rbf_func = "gaussian")
plot_epsilon_res(coexp)
get_epsilon_res(coexp)

coexp <- cor_module_coremo_clustering(coexp,
           coremo_params = params_coremo(epsilon = 2, rbf_func = "gaussian"))
coexp <- cor_module_coremo_stability(coexp)
coexp <- cor_module_coremo_cor_sign(coexp)
coexp <- cor_module_coremo_eigengene(coexp)

modules <- get_outputs(coexp)$final_modules
```

The epsilon trap: the scale-free R² criterion tends to pick a **high** epsilon,
and epsilon controls how hard the RBF suppresses weak correlations. A high value
throws genes out of modules entirely. Don't take the R² optimum on faith, scan a
few values and count how many features survive into a module.

### Graph / Leiden

```r
coexp <- cor_module_processing(coexp, cor_method = "spearman")
coexp <- cor_module_check_epsilon(coexp, rbf_func = "bump")

graph_params <- params_cor_graph(epsilon = 2)
coexp <- cor_module_graph_check_res(coexp,
           resolution_params = params_graph_resolution(),
           graph_params = graph_params)
plot_resolution_res(coexp)

coexp <- cor_module_graph_final_modules(coexp, min_size = 10L, max_size = 500L,
           subclustering = TRUE, .graph_params = graph_params)

modules <- get_modules(get_results(coexp))
```

### ICA

```r
coexp <- ica_processing(coexp, fast_svd = TRUE)

coexp <- ica_evaluate_comp(coexp,
           ica_type = "logcosh",
           ncomp_params = params_ica_ncomp(max_no_comp = 30L, steps = 2L),
           iter_params  = params_ica_randomisation(cross_validate = FALSE,
                                                   random_init = 20L))
plot_ica_ncomp_params(coexp)
coexp <- ica_optimal_ncomp(coexp)
coexp <- ica_stabilised_results(coexp)

res      <- get_results(coexp)
loadings <- get_factors(res, which = "gene_loadings")
modules  <- get_modules(res)
get_ica_stability_res(coexp)
```

`ica_optimal_ncomp()` fits a loess to the combined score against component
count and takes the inflection point. Give it a dense enough grid, on a handful
of points it gives up and warns. ICA modules **overlap**, so a gene can appear
in several.

### NMF (HALS)

Needs non-negative input, so build from CPM without the log transform:

```r
coexp <- BulkCoExp(raw_data = data_raw, meta_data = meta_data)
coexp <- preprocess_bulk_coexp(coexp, hvg = 0.5)

coexp <- nmf_bulk(coexp, k = 4L, preprocessing = "none",
                  nmf_hals_params = params_nmf_hals(max_iter = 300L, tol = 1e-4),
                  seed = 42L)

get_nmf_modules(coexp)          # long format, gene to module
get_nmf_gene_loadings(coexp)    # features x k
get_nmf_sample_activity(coexp)  # samples x k
```

`stabilised_nmf_bulk()` runs multiple random restarts and returns the one with
the lowest reconstruction loss. `get_nmf_stability()` reports the spread. Single
run with NNDSVD init is deterministic.

### DGRDL

```r
coexp <- BulkCoExp(raw_data = data_log, meta_data = meta_data) |>
  preprocess_bulk_coexp(hvg = 0.5, scaling = TRUE, scaling_type = "normal")

coexp <- dgrdl_grid_search(coexp, neighbours_vec = c(5L, 10L), ...)
get_grid_search_res(coexp)

coexp <- dgrdl_result(coexp, dgrdl_params = params_dgrdl())
get_modules(get_results(coexp))
```

Note the scaling in preprocessing, DGRDL wants it.

## Differential correlation

```r
coexp <- diffcor_module_processing(coexp, ...)
get_diffcor_graph(coexp)
```

Modules that change their correlation structure between conditions, rather than
their mean expression.

## Contrastive PCA

Finds structure in a target dataset that is not present in a background one.

```r
obj <- contrastive_pca_processing(target, background)
c_pca_plot_alphas(obj)                 # sweep the contrast strength
obj <- contrastive_pca(obj, alpha = ...)

get_c_pca_factors(obj)
get_c_pca_loadings(obj)
```

`alpha` is the knob: 0 is plain PCA on the target, large values maximally
suppress whatever the background shares. Look at the sweep before committing.

## Synthetic data for testing

`params_synthetic_bulk_rnaseq()` plus the generators in
`R/data_synthetic_bulk.R` give you count matrices with planted modules, which is
how you sanity-check a module-detection setup before pointing it at real data.
