---
name: bixverse
description: How to USE the bixverse R package for bioinformatics and computational biology analysis. Covers installation, the single cell suite (SingleCells, MetaCells, SingleCellsMultiModal, SingleCellsSubset), gene set enrichment and GSEA, bulk differential expression and co-expression modules, network diffusion, and ontology semantic similarity. Use this whenever a task involves bixverse, or any of its classes and generics (SingleCells, BulkDge, BulkCoExp, NetworkDiffusions, OntologySim, find_hvg_sc, calculate_pca_sc, gse_hypergeometric, params_* bundles), or when a user asks for single cell QC, clustering, marker genes, doublet detection, batch correction, metacells, GO enrichment, GSEA, pathway activity scoring, or co-expression module detection in R. Apply even when the package is not named but the user is working in an R session where bixverse is loaded.
---

# Using bixverse

bixverse wraps bioinformatics workflows around a Rust core via extendr. It's an
R package, but the heavy numerics happen in Rust, so it handles data sizes that
would kill a normal R session. A million cells on a laptop is the design target,
not a stretch goal.

The API is large: roughly 515 user-facing functions across single cell, gene set
enrichment, bulk RNAseq, graphs and ontologies. Do not guess at signatures.
Check `references/api-index.md` for whether something exists, then `?fn` for how
to call it.

## The three things to know before writing any code

**1. Tunables come in `params_*()` bundles.** Rather than functions with forty
arguments, bixverse passes one validated list:

```r
find_neighbours_sc(obj, neighbours_params = params_sc_neighbours(
  knn = list(k = 20L, knn_method = "hnsw")
))
```

Each bundle has a matching checkmate assertion, so a malformed one fails at the
call site with a readable message. There are 69 of them and they're all in the
API index.

**2. Analysis objects are S7 and the chain returns the object.** The shape is
`Constructor()` then `preprocess_*()` then verb methods then `get_results()`.
Everything pipes with `|>`. Two exceptions in single cell mutate the DuckDB in
place instead, see `references/single-cell.md`.

**3. Print the object.** Every class has a `print` method that tells you what
state it's in, what's been computed and what's gone stale. That's far cheaper
than guessing, and much cheaper than recomputing a PCA you didn't need to.

## Where to look

| Task | Read |
|---|---|
| Installing, Rust toolchain, missing dependencies, sister packages | `references/install.md` |
| `params_*()`, S7 chains, naming, getters, what mutates what | `references/conventions.md` |
| Single cell: load, QC, doublets, HVG, PCA, batch correction, clustering, embeddings, markers | `references/single-cell.md` |
| Single cell downstream: gene set scoring, SCENIC, topic models, differential expression, Hotspot, trajectory, miloR, metacells, CITE-seq, reference mapping | `references/single-cell-analysis.md` |
| Hypergeometric tests, GO elimination, GSEA, GSVA, ssGSEA, singscore | `references/enrichment.md` |
| Bulk RNAseq: differential expression, co-expression modules, ICA, NMF, contrastive PCA | `references/bulk.md` |
| Network diffusion, RBH graphs, similarity network fusion, ontology semantic similarity | `references/graphs-ontology.md` |
| Does function X exist? What's it called? | `references/api-index.md` |

## Smoke test

Runs on synthetic data, no downloads, finishes in seconds. Use it to confirm an
install works before debugging anything else.

```r
library(bixverse)

dir_sc <- file.path(tempdir(), "bixverse_smoke")
dir.create(dir_sc, showWarnings = FALSE)

syn <- generate_single_cell_test_data()

obj <- load_r_data(
  object = SingleCells(dir_data = dir_sc),
  counts = syn$counts,
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

obj  # print it, you should see the cache populated

unlink(dir_sc, recursive = TRUE, force = TRUE)
```

## Traps that bite everywhere

- **This is not Seurat or Scanpy.** Don't translate idioms across. Where an
  equivalent exists the name differs and the semantics usually differ too.
- **Single cell normalisation happens at load time**, driven by `target_size` in
  `params_sc_min_quality()`. There is no `normalise_sc()`. Looking for one and
  failing to find it does not mean you should write one.
- **A missing prerequisite warns and returns the object unchanged.** It does not
  error. `calculate_pca_sc()` without an HVG selection is a silent no-op. Check
  the object after each step rather than assuming a chain ran.
- **`rs_*` functions are raw extendr bindings with no input validation.** There
  is an R wrapper for each one that does the checking. Use that. The `rs_*`
  layer exists for people building on top of bixverse.
- **Indices handed to and from Rust are 0-based.** Getters like
  `get_cells_to_keep()` and `get_hvg()` return 0-based. Most setters take
  1-based or character names. Details in `references/single-cell.md`, and get
  this wrong and you silently analyse the wrong cells.
- **Old snake_case constructors are deprecated** since 0.3.0: `single_cell_exp`,
  `bulk_coexp`, `bulk_dge`, `network_diffusions`, `rbh_graph`, `snf`,
  `gene_ontology_data`, `ontology`, `meta_cells`. Use the UpperCamelCase S7 ones
  (`SingleCells`, `BulkCoExp`, ...). Older blog posts and forum answers still
  show the deprecated names.
- **Plotting mostly lives in `bixverse.plots`**, a separate package. bixverse
  itself gives you `extract_*_data()` functions that return the plot-ready
  data.table. If a plotting function isn't found, that's the missing dependency.
- **Windows is supported**, and is in CI. The MAX_PATH and cross-ABI handling
  for the HDF5 build lives in `tools/config.R`. Older answers saying otherwise
  predate that work.
