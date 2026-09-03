# bixverse conventions

Learn these once and most of the 475-function surface becomes predictable.

## Parameter bundles

Anything with more than a handful of tunables takes a single validated list
built by a `params_*()` constructor. There are 69 of them.

Naming:

- `params_<method>()` for a method's own tuning: `params_coremo()`,
  `params_scrublet()`, `params_snf()`, `params_gsea()`
- `params_sc_<method>()` when the method is single cell specific:
  `params_sc_harmony()`, `params_sc_hotspot()`, `params_sc_pca()`
- `params_*_defaults()` for zero-argument bundles that nest inside a bigger
  one: `params_knn_defaults()`, `params_hvg_defaults()`,
  `params_norm_doublets_defaults()`, `params_scenic_random_forest_defaults()`

Each returns a plain named list and gets checked by a matching `assert*()` from
the checkmate extensions, so a malformed bundle errors at the call site with a
message naming the offending element. Guessing an element name is safe in the
sense that it fails loudly rather than silently.

### The nesting quirk

Nested bundles are merged and then **flattened**. `params_sc_neighbours()` takes
a `knn` list, `modifyList`s it over `params_knn_defaults()`, and returns one
flat list:

```r
p <- params_sc_neighbours(knn = list(k = 20L, knn_method = "hnsw"))
p$k            # 20  <- flattened to the top level
p$knn$k        # NULL <- there is no nested $knn
```

So you pass `knn = list(...)` but read `p$k`. Same pattern in
`params_sc_fast_cluster()` and `params_scrublet()`.

A documented trap worth repeating: `pruning` in `params_sc_neighbours()` is a
share of the neighbourhood, so it has to track `k`. The default `1/12` is tuned
for the default `k = 15L`. Raise `k` without raising `pruning` and you prune
harder, which fails quietly: you still get a clustering, but under-connected
cells drop out as singleton communities and show up later as one-cell clusters
with inflated `run_paga_sc()` connectivities.

### Bundle or plain arguments?

Not uniform. These take plain arguments, no bundle:

`find_clusters_sc()`, `find_markers_sc()`, `find_all_markers_sc()`,
`umap_sc()` / `tsne_sc()` / `phate_sc()` (they take `manifoldsR::params_umap()`
and friends instead), `run_cell_qc()`, `set_cells_to_keep()`,
`gene_set_proportions_sc()`, `top_genes_perc_sc()`, `get_pseudobulked_sc()`,
`save_sc_exp_to_disk()`, `load_existing()`.

Nearly everything else takes one. When unsure, read `?fn`.

## Class chains

Analysis classes are S7 and inherit from `BixverseBaseClass`. The shape:

```r
obj <- ClassConstructor(data, ...)     # UpperCamelCase
obj <- preprocess_*(obj, ...)          # optional, domain dependent
obj <- verb_method(obj, params_*())    # repeat, each returns the object
res <- get_results(obj)                # the answer
```

Every step returns the object, so it all pipes:

```r
res <- BulkCoExp(raw_data, meta_data) |>
  preprocess_bulk_coexp(hvg = 0.5) |>
  cor_module_processing() |>
  get_results()
```

The universal getters are `get_results()`, `get_params()`, `get_metadata()` and
`get_outputs()`. Domain-specific getters exist alongside them
(`get_modules()`, `get_factors()`, `get_dge_list()`, `get_sim_matrix()`, ...).

## Verb vocabulary

| Prefix | Contract |
|---|---|
| `get_` | reads, never mutates |
| `set_` | mutates one slot |
| `calculate_` / `calc_` | computes and stores back on the object |
| `find_` | discovery: HVGs, neighbours, clusters, markers |
| `run_` | runs a named algorithm end to end |
| `generate_` | builds a derived object (kNN, RBH graph, metacells) |
| `add_` | appends data, a modality, or metadata |
| `plot_` | returns a ggplot |
| `extract_*_data` | returns a plot-ready data.table for `bixverse.plots` |
| `rs_` | raw extendr binding, no validation, not for normal use |

`calculate_` and `calc_` mean the same thing. The inconsistency is historical.

## Printing beats guessing

Every class and every result object has a `print` method. Printing a
`SingleCells` shows dimensions, the filter state, what's in the cache and what's
gone stale. Do that before deciding a step needs re-running.

```r
sc_object          # state of the object
qc                 # a CellQc summary
get_params(obj)    # exactly what was used
```

## Deprecated constructors

The snake_case constructors were replaced by S7 UpperCamelCase ones in 0.3.0.
They still work and emit a lifecycle warning, but older material and forum
answers use them, so recognise them:

| Deprecated | Use |
|---|---|
| `single_cell_exp()` | `SingleCells()` |
| `meta_cells()` | `MetaCells()` |
| `bulk_coexp()` | `BulkCoExp()` |
| `bulk_dge()` | `BulkDge()` |
| `network_diffusions()` | `NetworkDiffusions()` |
| `rbh_graph()` | `RbhGraph()` |
| `snf()` | `SimilarityNetworkFusion()` |
| `gene_ontology_data()` | `GeneOntologyElim()` |
| `ontology()` | `OntologySim()` |

## data.table everywhere

Results come back as `data.table`, not tibbles or data.frames. Chain with
`[i, j, by]` rather than converting. The package imports data.table wholesale,
so the syntax is available once bixverse is attached.
