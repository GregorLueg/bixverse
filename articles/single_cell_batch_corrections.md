# Batch correction with bixverse

## Intro

Batch effects are systematic technical differences between experimental
runs that can obscure genuine biological signal in single cell data.
When datasets from multiple experiments, technologies, or laboratories
are combined, cells tend to group by their source rather than by cell
type. Batch correction methods aim to remove these technical artefacts
whilst preserving biological variation.

`bixverse` provides four batch correction approaches, each with
different trade-offs:

- **fastMNN** ([Haghverdi, et al.,
  2018](https://doi.org/10.1038/nbt.4091)) identifies mutual nearest
  neighbours across batches: pairs of cells from different batches that
  are each other’s closest match. Correction vectors are computed from
  these pairs and applied to align the embeddings. It is particularly
  effective when cell type composition differs across batches.

- **Harmony** ([Korsunsky, et al.,
  2019](https://doi.org/10.1038/s41592-019-0619-0)) operates on PCA
  embeddings using iterative soft clustering with diversity penalties.
  It assigns cells to clusters, estimates batch effects per cluster via
  ridge regression, and corrects the embedding. It tends to be fast and
  works well across a range of scenarios. `bixverse` defaults to version
  2 ([Patikas, et al.,
  2026](https://www.biorxiv.org/content/10.64898/2026.03.16.711825v1)),
  which scales better to large data sets and changes the objective so it
  doesn’t overcorrect when the batches hold different cell types. The
  original is still there for reproducing old results.

- **BBKNN** ([Polanski, et al.,
  2020](https://doi.org/10.1093/bioinformatics/btz625)) takes a
  fundamentally different approach. Rather than correcting an embedding,
  it constructs a batch-balanced k-nearest neighbour graph by building
  separate neighbour indices per batch. The corrected graph can then be
  used directly for clustering and visualisation.

- **Seurat anchors** ([Stuart, et al.,
  2019](https://doi.org/10.1016/j.cell.2019.05.031)) find pairs of cells
  that are mutual nearest neighbours in a space shared by two batches,
  score them by how much of their neighbourhood they share, and use them
  to pull one batch onto the other. Two flavours are available. CCA
  builds the shared space from canonical correlations and corrects
  hardest, which suits batches that don’t share their cell types. rPCA
  projects each batch into the other’s PCA basis instead: cheaper,
  gentler, and the better default when the batches broadly match. Both
  need batch-aware HVGs.

A correction has two jobs: mix the batches and keep the biology intact.
Going too hard on the first ruins the second, so `bixverse` scores both,
following the scIB benchmark ([Luecken, et al.,
2022](https://doi.org/10.1038/s41592-021-01336-8)).

Batch mixing:

- **kBET**: tests whether neighbourhood batch proportions match the
  global ones.
- **Batch ASW**: silhouette width on batch labels in the embedding.
- **iLISI**: effective number of batches per neighbourhood. Works on any
  kNN graph.
- **PCR**: how much of the embedding variance batch explains, before
  versus after correction.

Bio conservation, which needs cell type labels:

- **cLISI**: effective number of cell types per neighbourhood. You want
  one.
- **Cell type ASW**: silhouette width on cell type labels.
- **Graph connectivity**: whether each cell type stays one connected
  piece of the kNN graph.

[`calculate_integration_metrics_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_integration_metrics_sc.md)
runs the lot and returns one row, every column rescaled so that higher
is better.

``` r

library(bixverse)
library(bixverse.plots)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
library(ggplot2)
library(magrittr)
```

## Preparing the data

### Loading the batches

We use two PBMC datasets (pbmc3k and pbmc4k) as a standard benchmark for
batch correction. These come from the same tissue but differ in
sequencing depth and cell counts, introducing a clear batch effect.

Code

``` r

dir_data <- download_pbmc_batches()

tempdir_batch_cor <- file.path(tempdir(), "batch_cor_bixverse")
dir.create(tempdir_batch_cor, showWarnings = FALSE, recursive = TRUE)

h5ad_files <- list.files(dir_data)
h5ad_files <- h5ad_files[grepl(".h5ad", h5ad_files)]
h5ad_paths <- file.path(dir_data, h5ad_files)
names(h5ad_paths) <- c("pbmc3k", "pbmc4k")

h5_tasks <- prescan_h5ad_files(h5_paths = h5ad_paths)

sc_object <- SingleCells(dir_data = tempdir_batch_cor)

sc_object <- load_multi_h5ad(
  object = sc_object,
  prescan_result = h5_tasks,
  .verbose = TRUE
)
#>  Converting the cell-based data into the gene-based format.
#> Loading observation data from h5ad files into DuckDB.
#> Loading variable data into DuckDB.
```

### Quality control

Before batch correction, we filter out low-quality cells using
mitochondrial gene proportions and library complexity metrics.

Code

``` r

var <- get_sc_var(sc_object)

# let's get the gene symbols from one of the h5ad files - in multi file import
# additional columns are dropped from the vars
h5_metadata <- read_h5ad_metadata(h5ad_paths[[1]])

var <- merge(
  var,
  h5_metadata$var[, c("ENSEMBL_ID", "Symbol_TENx")],
  by.x = "gene_id",
  by.y = "ENSEMBL_ID"
)

setnames(var, old = "Symbol_TENx", new = "gene_symbol", skip_absent = TRUE)

gs_of_interest <- list(
  MT = var[grepl("^MT-", gene_symbol), gene_id],
  Ribo = var[grepl("^RPS|^RPL", gene_symbol), gene_id]
)

sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  streaming = FALSE,
  .verbose = TRUE
)

qc_df <- sc_object[[c("cell_id", "lib_size", "nnz", "MT")]]

metrics <- list(
  log10_lib_size = log10(qc_df$lib_size),
  log10_nnz = log10(qc_df$nnz),
  MT = qc_df$MT
)
directions <- c(
  log10_lib_size = "twosided",
  log10_nnz = "twosided",
  MT = "above"
)

qc <- run_cell_qc(
  metrics = metrics,
  cells_to_keep = get_cells_to_keep(sc_object),
  directions = directions,
  threshold = 3
)

# Set the cells to keep
sc_object[["outlier"]] <- qc$combined
cells_to_keep <- qc_df[!qc$combined, cell_id]
sc_object <- set_cells_to_keep(sc_object, cells_to_keep)
```

### Pre-processing

Standard pre-processing: HVG selection, PCA, and neighbour computation.

Code

``` r

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = 2000L,
  .verbose = TRUE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 30L
)
#> Using sparse SVD solving on scaled data on 2000 HVG.

sc_object <- find_neighbours_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours()
)
#> 
#> Generating the snn graph.
```

### Cell type labels

The bio conservation metrics need cell types. We cluster the uncorrected
data and annotate the clusters with scType ([Ianevski, et al.,
2022](https://doi.org/10.1038/s41467-022-28803-w)) on a handful of
canonical PBMC markers, same as in the PBMC vignette. The labels are
deliberately coarse: T cells stay lumped, since splitting CD4 from CD8
on uncorrected data is asking for trouble.

One caveat: the labels come from uncorrected data, so they carry
whatever the batch effect did to the clusters. Good enough to spot a
method that smears cell types together, not a ground truth.

Code

``` r

sc_object <- find_clusters_sc(sc_object, res = 1, name = "leiden_uncorrected")
#> Running leiden clustering at resolution 1.

cell_markers <- c(
  CD3D = "T cells",
  CD3E = "T cells",
  IL7R = "T cells",
  CD8A = "T cells",
  MS4A1 = "B cells",
  CD79A = "B cells",
  CD14 = "Monocytes",
  LYZ = "Monocytes",
  FCGR3A = "Monocytes",
  GNLY = "NK",
  NKG7 = "NK",
  FCER1A = "DC",
  CD1C = "DC",
  PPBP = "Platelets"
)

cell_markers_dt <- data.table(
  gene_symbol = names(cell_markers),
  cell_type = unname(cell_markers)
) %>%
  .[, gene_id := var$gene_id[match(gene_symbol, var$gene_symbol)]] %>%
  .[!is.na(gene_id)]

marker_list <- prepare_cell_markers(sc_object, cell_markers_dt)
sctype_scores <- calc_sc_type_scores(
  object = sc_object,
  cell_marker_list = marker_list
)

cell_type_anno <- score_clusters(
  sctype_scores,
  sc_object[[]][["leiden_uncorrected"]]
)

obs <- get_sc_obs(sc_object, filtered = TRUE)[, .(
  cell_idx,
  leiden_uncorrected
)] %>%
  .[,
    cell_type := cell_type_anno$cell_type[match(
      leiden_uncorrected,
      cell_type_anno$cluster_id
    )]
  ]

sc_object[["cell_type"]] <- obs$cell_type

table(obs$cell_type)
#> 
#>   B cells        DC Monocytes        NK Platelets   T cells 
#>       933       163      1421       889        28      3192
```

## Comparing the different methods

Let’s check out the different methods and how they behave.

### Uncorrected data

Before applying any correction, we can visualise the batch effect and
compute baseline metrics.

``` r

sc_object <- tsne_sc(object = sc_object, approx_type = "bh_qd")
#> Running t-SNE.
```

tSNE plot (tSNE’s usually look prettier… It’s anyway not a good way to
assess batch effect correction, so at least let’s make it pleasant on
the eye. We use the quick-and-dirty version which is fast).

``` r

embedding_plot_sc(
  sc_object,
  embedding = "tsne",
  colour_by = "exp_id",
  label_by = "exp_id",
  discrete = TRUE
) +
  labs(
    title = "Before batch correction",
    colour = "Batch:"
  )
```

![](single_cell_batch_corrections_files/figure-html/uncorrected%20plot-1.png)

Each metric has its own function with its own print method. Here’s kBET
on its own:

``` r

calculate_kbet_sc(sc_object, batch_column = "exp_id")
#> kBET Scores
#>   Cells: 6626 | Batches: 2 | Threshold: 0.050
#>   Rejection rate:      0.9870 (6540 / 6626)
#>   Mean Chi-Square:     12.7892 (expected under H0: 1)
#>   Median Chi-Square:   8.1742
```

In practice, you’d want all of them at once.
[`calculate_integration_metrics_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_integration_metrics_sc.md)
gives one row per call, and we’ll stack one row per method into
`metrics_dt`:

``` r

metrics_dt <- calculate_integration_metrics_sc(
  sc_object,
  batch_column = "exp_id",
  cell_type_column = "cell_type",
  .verbose = FALSE
) %>%
  .[, method := "Uncorrected"]

metrics_dt[]
#>    embedding kbet_accept batch_asw ilisi pcr_comparison clisi cell_type_asw
#>       <char>       <num>     <num> <num>          <num> <num>         <num>
#> 1:       pca  0.01297917 0.9063555     0             NA     1     0.6508355
#>    graph_connectivity      method
#>                 <num>      <char>
#> 1:          0.9501534 Uncorrected
```

The high kBET rejection rate (so a low `kbet_accept`) and a low `ilisi`
confirm a substantial batch effect, and the tSNE shows clear separation
by batch. `pcr_comparison` is `NA` here: it compares a corrected
embedding against this PCA, so the PCA has nothing to compare against.

### fastMNN

fastMNN works on batch-aware HVGs and produces a corrected embedding. It
identifies mutual nearest neighbours across batches to compute
correction vectors. For it to work best, we need the batch-aware HVGs.
We will calculate these, provide them to fastMNN (which regenerates the
PCA embedding based on the batch-aware HVG) and run then the actual
algorithm

``` r

batch_aware_hvg <- find_hvg_batch_aware_sc(
  object = sc_object,
  batch_column = "exp_id"
)

sc_object <- fast_mnn_sc(
  object = sc_object,
  batch_hvg_genes = batch_aware_hvg$hvg_gene_idx,
  batch_column = "exp_id"
)

sc_object <- find_neighbours_sc(
  object = sc_object,
  embd_to_use = "mnn",
  neighbours_params = params_sc_neighbours()
)
#> 
#> Generating the snn graph.
```

``` r

metrics_dt <- rbind(
  metrics_dt,
  calculate_integration_metrics_sc(
    sc_object,
    batch_column = "exp_id",
    cell_type_column = "cell_type",
    embd_to_use = "mnn",
    .verbose = FALSE
  ) %>%
    .[, method := "fastMNN"]
)

metrics_dt[method == "fastMNN"]
#>    embedding kbet_accept batch_asw     ilisi pcr_comparison clisi cell_type_asw
#>       <char>       <num>     <num>     <num>          <num> <num>         <num>
#> 1:       mnn   0.6845759 0.9736488 0.4705881      0.9600038     1     0.6492018
#>    graph_connectivity  method
#>                 <num>  <char>
#> 1:          0.9489644 fastMNN
```

We can see clear improvements on the batch side: kBET acceptance, batch
ASW and iLISI all go up.

``` r

sc_object <- tsne_sc(
  object = sc_object,
  slot_name = "tsne_mnn",
  use_knn = TRUE,
  approx_type = "bh_qd"
)
#> Running t-SNE.
```

``` r

embedding_plot_sc(
  sc_object,
  embedding = "tsne_mnn",
  colour_by = "exp_id",
  label_by = "exp_id",
  discrete = TRUE
) +
  labs(
    title = "fastMNN batch correction",
    colour = "Batch:"
  )
```

![](single_cell_batch_corrections_files/figure-html/fastmnn%20plot-1.png)

Also, visually the batches get mixed now.

### Harmony

Harmony operates on the PCA embedding directly. Leave `k` as `NULL` and
the number of clusters is set to one per 30 cells, capped at 100. Both
versions are supported and share one numerical core in Rust.
[`harmony_sc()`](https://gregorlueg.github.io/bixverse/reference/harmony_sc.md)
runs v2 and stores `"harmony"`,
[`harmony_v1_sc()`](https://gregorlueg.github.io/bixverse/reference/harmony_v1_sc.md)
runs the original and stores `"harmony_v1"`. The intercept is never
penalised in either.

#### Version 2 (default)

[Patikas, et al.,
2026](https://www.biorxiv.org/content/10.64898/2026.03.16.711825v1)
changed four things compared to v1. The last three have knobs in
[`params_sc_harmony_v2()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony_v2.md):

- Stabilised diversity penalty, which makes the objective
  scale-invariant.
- Batch pruning in the ridge regression. A batch that puts less than
  `batch_proportion_cutoff` of its cells into a cluster is dropped from
  that cluster’s regression. A cluster left with a single batch isn’t
  corrected, so a cell type unique to one batch doesn’t get pulled
  towards the others.
- Dynamic lambda. With `use_dynamic_lambda = TRUE` (the default) the
  ridge penalty is set per cluster as `alpha` times the expected counts
  and the fixed `lambda` is ignored.
- Theta scaling by batch size via `tau`. Off by default (`tau = 0`).

The paper’s arrowhead matrix inversion for a single batch covariate is
used by both versions here, so v1 gets that speed-up too. v2 also runs
only 4 k-means iterations per Harmony round versus 20 in v1.

``` r

sc_object <- harmony_sc(
  object = sc_object,
  batch_column = "exp_id",
  harmony_params = params_sc_harmony_v2()
)
#>  Auto-determined number of Harmony clusters: 100

sc_object <- find_neighbours_sc(
  object = sc_object,
  embd_to_use = "harmony",
  neighbours_params = params_sc_neighbours()
)
#> 
#> Generating the snn graph.

# for the cluster comparison in the BBKNN section
sc_object <- find_clusters_sc(sc_object, res = 1, name = "leiden_harmony")
#> Running leiden clustering at resolution 1.
```

Let’s calculate the batch correction-related scores

``` r

metrics_dt <- rbind(
  metrics_dt,
  calculate_integration_metrics_sc(
    sc_object,
    batch_column = "exp_id",
    cell_type_column = "cell_type",
    embd_to_use = "harmony",
    .verbose = FALSE
  ) %>%
    .[, method := "Harmony v2"]
)

metrics_dt[method == "Harmony v2"]
#>    embedding kbet_accept batch_asw     ilisi pcr_comparison clisi cell_type_asw
#>       <char>       <num>     <num>     <num>          <num> <num>         <num>
#> 1:   harmony   0.7982191  0.977846 0.6423357      0.9493217     1     0.6677088
#>    graph_connectivity     method
#>                 <num>     <char>
#> 1:          0.9499137 Harmony v2
```

Also here, we observe improvements across the board.

``` r

sc_object <- tsne_sc(
  object = sc_object,
  slot_name = "tsne_harm",
  use_knn = TRUE,
  approx_type = "bh_qd"
)
#> Running t-SNE.
```

``` r

embedding_plot_sc(
  sc_object,
  embedding = "tsne_harm",
  colour_by = "exp_id",
  label_by = "exp_id",
  discrete = TRUE
) +
  labs(
    title = "Harmony (version 2) batch correction",
    colour = "Batch:"
  )
```

![](single_cell_batch_corrections_files/figure-html/harmony%20plot-1.png)

Got a GPU? `harmony_v2_gpu_sc()` in
[bixverse.gpu](https://gregorlueg.github.io/bixverse.gpu/articles/gpu_single_cell.html)
runs the same v2 algorithm and stores the result as `"harmony_gpu"`. It
takes a single batch column only.

#### Version 1

Need to reproduce an old result? The original lives on as
[`harmony_v1_sc()`](https://gregorlueg.github.io/bixverse/reference/harmony_v1_sc.md),
with its parameters in
[`params_sc_harmony()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony.md).
Otherwise, stick with v2.

``` r

sc_object <- harmony_v1_sc(
  object = sc_object,
  batch_column = "exp_id",
  harmony_params = params_sc_harmony()
)
#>  Auto-determined number of Harmony clusters: 100

sc_object <- find_neighbours_sc(
  object = sc_object,
  embd_to_use = "harmony_v1",
  neighbours_params = params_sc_neighbours()
)
#> 
#> Generating the snn graph.
```

``` r

metrics_dt <- rbind(
  metrics_dt,
  calculate_integration_metrics_sc(
    sc_object,
    batch_column = "exp_id",
    cell_type_column = "cell_type",
    embd_to_use = "harmony_v1",
    .verbose = FALSE
  ) %>%
    .[, method := "Harmony v1"]
)

metrics_dt[method == "Harmony v1"]
#>     embedding kbet_accept batch_asw     ilisi pcr_comparison clisi
#>        <char>       <num>     <num>     <num>          <num> <num>
#> 1: harmony_v1   0.8827347 0.9754045 0.6423357      0.9817465     1
#>    cell_type_asw graph_connectivity     method
#>            <num>              <num>     <char>
#> 1:     0.6671693          0.9499659 Harmony v1
```

``` r

sc_object <- tsne_sc(
  object = sc_object,
  slot_name = "tsne_harm_v1",
  use_knn = TRUE,
  approx_type = "bh_qd"
)
#> Running t-SNE.
```

``` r

embedding_plot_sc(
  sc_object,
  embedding = "tsne_harm_v1",
  colour_by = "exp_id",
  label_by = "exp_id",
  discrete = TRUE
) +
  labs(
    title = "Harmony (version 1) batch correction",
    colour = "Batch:"
  )
```

![](single_cell_batch_corrections_files/figure-html/harmony%20(v1)%20plot-1.png)

### Seurat CCA

Seurat’s anchor methods take a different route. For each pair of batches
they find cells that are mutual nearest neighbours in a shared space,
treat those pairs as anchors, score them by how much of their
neighbourhood they share, and then apply a kernel-weighted correction
that pulls the query batch onto the reference. Batches get merged in the
order of their pairwise anchor counts.

For CCA the shared space is a canonical correlation embedding computed
per batch pair. Anchors also get filtered in gene space before scoring,
which is CCA-only in Seurat.

A few things differ from Seurat’s reference implementation:

- The per-gene `ScaleData` step is skipped. We work directly from
  per-cell standardised log-normalised HVG expression.
- `M = X1^T @ X2` is never materialised. The canonical correlations come
  from a matrix-free randomised SVD, so memory stays at
  `O((n1 + n2) * num_cc)` rather than `n1 x n2`. On big batch pairs that
  is the difference between running and not running.
- The correction is applied to the embedding (dims x cells), not to full
  log-expression.
- The effective CC rank is `max(num_cc, dims)`, not `num_cc`.

Anchor structure comes out close to identical. Like fastMNN, CCA wants
batch-aware HVGs, so we reuse the ones computed above.

``` r

sc_object <- seurat_cca_sc(
  object = sc_object,
  batch_hvg_genes = batch_aware_hvg$hvg_gene_idx,
  batch_column = "exp_id"
)

sc_object <- find_neighbours_sc(
  object = sc_object,
  embd_to_use = "cca",
  neighbours_params = params_sc_neighbours()
)
#> 
#> Generating the snn graph.
```

``` r

metrics_dt <- rbind(
  metrics_dt,
  calculate_integration_metrics_sc(
    sc_object,
    batch_column = "exp_id",
    cell_type_column = "cell_type",
    embd_to_use = "cca",
    .verbose = FALSE
  ) %>%
    .[, method := "Seurat CCA"]
)

metrics_dt[method == "Seurat CCA"]
#>    embedding kbet_accept batch_asw     ilisi pcr_comparison clisi cell_type_asw
#>       <char>       <num>     <num>     <num>          <num> <num>         <num>
#> 1:       cca     0.77377 0.9569454 0.4705881      0.9784447     1     0.6653001
#>    graph_connectivity     method
#>                 <num>     <char>
#> 1:          0.9499659 Seurat CCA
```

``` r

sc_object <- tsne_sc(
  object = sc_object,
  slot_name = "tsne_cca",
  use_knn = TRUE,
  approx_type = "bh_qd"
)
#> Running t-SNE.
```

``` r

embedding_plot_sc(
  sc_object,
  embedding = "tsne_cca",
  colour_by = "exp_id",
  label_by = "exp_id",
  discrete = TRUE
) +
  labs(
    title = "Seurat CCA batch correction",
    colour = "Batch:"
  )
```

![](single_cell_batch_corrections_files/figure-html/seurat%20cca%20plot-1.png)

### Seurat rPCA

rPCA runs the same anchor pipeline on a cheaper shared space. Instead of
computing canonical correlations, each batch keeps its own PCA basis and
the other batch’s expression is projected into it. Cross-batch neighbour
queries then happen in these projected bases.

The trade-off is the usual one: rPCA is faster and corrects less
aggressively, which makes it the safer choice when your batches share
most of their cell types. CCA is the one to reach for when they don’t.
As in Seurat, no gene-space anchor filter runs here, that step is
CCA-only.

``` r

sc_object <- seurat_rpca_sc(
  object = sc_object,
  batch_hvg_genes = batch_aware_hvg$hvg_gene_idx,
  batch_column = "exp_id"
)

sc_object <- find_neighbours_sc(
  object = sc_object,
  embd_to_use = "rpca",
  neighbours_params = params_sc_neighbours()
)
#> 
#> Generating the snn graph.
```

``` r

metrics_dt <- rbind(
  metrics_dt,
  calculate_integration_metrics_sc(
    sc_object,
    batch_column = "exp_id",
    cell_type_column = "cell_type",
    embd_to_use = "rpca",
    .verbose = FALSE
  ) %>%
    .[, method := "Seurat rPCA"]
)

metrics_dt[method == "Seurat rPCA"]
#>    embedding kbet_accept batch_asw     ilisi pcr_comparison clisi cell_type_asw
#>       <char>       <num>     <num>     <num>          <num> <num>         <num>
#> 1:      rpca   0.8131603 0.9678163 0.6423357      0.9419804     1      0.663178
#>    graph_connectivity      method
#>                 <num>      <char>
#> 1:          0.9495387 Seurat rPCA
```

``` r

sc_object <- tsne_sc(
  object = sc_object,
  slot_name = "tsne_rpca",
  use_knn = TRUE,
  approx_type = "bh_qd"
)
#> Running t-SNE.
```

``` r

embedding_plot_sc(
  sc_object,
  embedding = "tsne_rpca",
  colour_by = "exp_id",
  label_by = "exp_id",
  discrete = TRUE
) +
  labs(
    title = "Seurat rPCA batch correction",
    colour = "Batch:"
  )
```

![](single_cell_batch_corrections_files/figure-html/seurat%20rpca%20plot-1.png)

One caveat worth knowing about both anchor methods: on data that already
mixes well they have nothing to gain, and the correction can nudge the
batch metrics slightly in the wrong direction while still sharpening the
cell type structure. Check your uncorrected metrics before reaching for
them.

### BBKNN

BBKNN doesn’t touch the embedding. It builds a batch-balanced neighbour
graph instead: `neighbours_within_batch` nearest neighbours from every
batch, turned into UMAP-style connectivities. Here that’s 10 per batch,
so 20 neighbours per cell feed the connectivities.
`no_neighbours_to_keep` only cuts the stored kNN, and the cut stays
balanced: each batch keeps its nearest `ceil(15 / 2) = 8`, so a cell
ends up with 8 from one batch and 7 from the other.

That makes BBKNN awkward to score. kBET and iLISI read the stored kNN,
which is balanced by construction, so they come out near perfect
whatever the data looks like. That’s the method grading its own
homework, not a sign of mixing. ASW and PCR need an embedding, and BBKNN
doesn’t have one, so we pass `embd_to_use = NULL` and those come back
`NA`. cLISI and graph connectivity are still fair, since they only look
at cell types.

``` r

sc_object <- bbknn_sc(
  object = sc_object,
  batch_column = "exp_id",
  no_neighbours_to_keep = 15L,
  bbknn_params = params_sc_bbknn(neighbours_within_batch = 10L)
)
#> Warning in `method(bbknn_sc, bixverse::SingleCells)`(object = <object>, : Prior
#> kNN matrix found. Will be overwritten.
#> Running BBKNN algorithm.
#> Generating graph based on BBKNN connectivities. Weights will be based on the connectivities and not shared nearest neighbour calculations.
```

``` r

metrics_dt <- rbind(
  metrics_dt,
  calculate_integration_metrics_sc(
    sc_object,
    batch_column = "exp_id",
    cell_type_column = "cell_type",
    embd_to_use = NULL,
    .verbose = FALSE
  ) %>%
    .[, method := "BBKNN"]
)

metrics_dt[method == "BBKNN"]
#>    embedding kbet_accept batch_asw     ilisi pcr_comparison clisi cell_type_asw
#>       <char>       <num>     <num>     <num>          <num> <num>         <num>
#> 1:      <NA>           1        NA 0.9911504             NA     1            NA
#>    graph_connectivity method
#>                 <num> <char>
#> 1:          0.9501534  BBKNN
```

So what does the graph actually look like? Counting edges is not enough:
the connectivities are weighted, and everything downstream (clustering,
t-SNE, UMAP) follows the weights. Below: the cross-batch share of the
stored kNN, of the graph’s edges and of its total weight, plus the
median distance from each cell to its nearest neighbour in its own batch
and in the other one.

``` r

batch <- get_sc_obs(sc_object, filtered = TRUE)$exp_id

# kNN indices are 1-based on the R side
knn_idx <- get_knn_mat(sc_object) + 1L
knn_dist <- get_knn_dist(sc_object)
other_batch <- matrix(batch[knn_idx] != batch, nrow = nrow(knn_idx))

graph <- as_igraph(get_connectivities(sc_object))
edge_ends <- igraph::ends(graph, igraph::E(graph), names = FALSE)
cross_edge <- batch[edge_ends[, 1]] != batch[edge_ends[, 2]]
weights <- igraph::E(graph)$weight

data.table(
  cross_batch_share_knn = mean(other_batch),
  cross_batch_share_edges = mean(cross_edge),
  cross_batch_share_weight = sum(weights[cross_edge]) / sum(weights),
  median_nearest_same_batch = median(
    matrixStats::rowMins(ifelse(other_batch, Inf, knn_dist))
  ),
  median_nearest_other_batch = median(
    matrixStats::rowMins(ifelse(other_batch, knn_dist, Inf))
  )
)
#>    cross_batch_share_knn cross_batch_share_edges cross_batch_share_weight
#>                    <num>                   <num>                    <num>
#> 1:             0.4668377               0.5072296               0.03727733
#>    median_nearest_same_batch median_nearest_other_batch
#>                        <num>                      <num>
#> 1:                  4.547937                   8.465939
```

Half the edges cross batches, but they carry only a few percent of the
weight. The nearest cell in the other batch sits about twice as far away
as the nearest one in the cell’s own batch, and the UMAP connectivities
decay with distance, so the cross-batch edges end up with next to no
weight. This isn’t an implementation detail: the reference `bbknn`
package behaves the same on this data. pbmc3k and pbmc4k were run on
different 10x chemistries (GemCode v1 against v2), and that gap is real.

The fair test of a graph method is clustering, so let’s run Leiden on
the BBKNN graph and compare its clusters against the uncorrected ones
and Harmony v2. A cluster where one batch holds more than 90% of the
cells counts as unmixed.

``` r

sc_object <- find_clusters_sc(sc_object, res = 1, name = "leiden_bbknn")
#> Running leiden clustering at resolution 1.

# share of a cluster's cells from one batch above which it counts as unmixed
single_batch_cutoff <- 0.9

batch_composition <- function(cluster_col) {
  obs <- get_sc_obs(sc_object, filtered = TRUE)[, .(
    cluster = get(cluster_col),
    exp_id
  )]
  per_cluster <- obs[, .N, .(cluster, exp_id)] %>%
    .[, share := N / sum(N), cluster] %>%
    .[, .(n_cells = sum(N), top_batch_share = max(share)), cluster]
  unmixed <- per_cluster[top_batch_share > single_batch_cutoff]
  data.table(
    clustering = cluster_col,
    n_clusters = nrow(per_cluster),
    n_unmixed_clusters = nrow(unmixed),
    share_cells_in_unmixed = sum(unmixed$n_cells) / nrow(obs)
  )
}

rbindlist(purrr::map(
  c("leiden_uncorrected", "leiden_harmony", "leiden_bbknn"),
  batch_composition
))
#>            clustering n_clusters n_unmixed_clusters share_cells_in_unmixed
#>                <char>      <int>              <int>                  <num>
#> 1: leiden_uncorrected         22                 17            0.971174162
#> 2:     leiden_harmony         15                  1            0.005131301
#> 3:       leiden_bbknn         18                 13            0.876546936
```

``` r

sc_object <- tsne_sc(object = sc_object, use_knn = TRUE, approx_type = "bh_qd")
#> Running t-SNE.
```

``` r

embedding_plot_sc(
  sc_object,
  embedding = "tsne",
  colour_by = "exp_id",
  label_by = "exp_id",
  discrete = TRUE
) +
  labs(
    title = "BBKNN batch correction",
    colour = "Batch:"
  )
```

![](single_cell_batch_corrections_files/figure-html/bbknn%20plot-1.png)

The clusters and the t-SNE tell the same story. Leiden on the BBKNN
graph still puts most cells into single-batch clusters, barely better
than doing nothing, while Harmony mixes nearly all of them. The t-SNE
follows the same weights and stays split by batch. On a batch effect
this strong, BBKNN on its own doesn’t correct much. It’s built for
milder effects, where the cross-batch neighbours sit close enough to
keep real weight.

## Choosing a method

Here’s everything side by side:

``` r

metrics_dt[, .(
  method,
  kbet_accept,
  batch_asw,
  ilisi,
  pcr_comparison,
  clisi,
  cell_type_asw,
  graph_connectivity
)]
#>         method kbet_accept batch_asw     ilisi pcr_comparison clisi
#>         <char>       <num>     <num>     <num>          <num> <num>
#> 1: Uncorrected  0.01297917 0.9063555 0.0000000             NA     1
#> 2:     fastMNN  0.68457591 0.9736488 0.4705881      0.9600038     1
#> 3:  Harmony v2  0.79821914 0.9778460 0.6423357      0.9493217     1
#> 4:  Harmony v1  0.88273468 0.9754045 0.6423357      0.9817465     1
#> 5:  Seurat CCA  0.77377000 0.9569454 0.4705881      0.9784447     1
#> 6: Seurat rPCA  0.81316028 0.9678163 0.6423357      0.9419804     1
#> 7:       BBKNN  1.00000000        NA 0.9911504             NA     1
#>    cell_type_asw graph_connectivity
#>            <num>              <num>
#> 1:     0.6508355          0.9501534
#> 2:     0.6492018          0.9489644
#> 3:     0.6677088          0.9499137
#> 4:     0.6671693          0.9499659
#> 5:     0.6653001          0.9499659
#> 6:     0.6631780          0.9495387
#> 7:            NA          0.9501534
```

``` r

bio_metrics <- c("clisi", "cell_type_asw", "graph_connectivity")

metrics_long <- melt(
  metrics_dt[, -"embedding"],
  id.vars = "method",
  variable.name = "metric",
  value.name = "score"
) %>%
  # BBKNN's kBET and iLISI read its own balanced kNN, so they are greyed out
  .[method == "BBKNN" & metric %in% c("kbet_accept", "ilisi"), score := NA] %>%
  .[, `:=`(
    method = factor(method, levels = rev(unique(metrics_dt$method))),
    type = fifelse(metric %in% bio_metrics, "Bio conservation", "Batch mixing"),
    label_col = fifelse(!is.na(score) & score > 0.6, "white", "black")
  )]

ggplot(metrics_long, aes(x = metric, y = method, fill = score)) +
  geom_tile(colour = "white") +
  geom_text(
    aes(
      label = fifelse(is.na(score), "", sprintf("%.2f", score)),
      colour = label_col
    ),
    size = 3
  ) +
  scale_colour_identity() +
  facet_grid(~type, scales = "free_x", space = "free_x") +
  scale_fill_distiller(
    palette = "Blues",
    direction = 1,
    na.value = "grey90",
    limits = c(0, 1),
    oob = scales::squish
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "Sc metrics")
```

![](single_cell_batch_corrections_files/figure-html/metrics%20plot-1.png)

Read the two panels together. A method that tops batch mixing but drops
on bio conservation has merged cell types along with the batches, which
is worse than doing nothing.

Two of these tell you little on this data set. cLISI sits at 1 for every
method: six coarse, well separated cell types keep their neighbourhoods
pure no matter what. Batch ASW barely moves either, since two PBMC
batches of the same tissue never separate much on silhouette to begin
with. Both earn their keep on finer labels and nastier batch effects. On
two PBMC runs, kBET, iLISI and PCR do the talking. BBKNN’s row is greyed
out on the batch side: its kBET and iLISI in the table above are the
artefact from its section, not a win.

There is no universally best batch correction method. Some practical
guidance:

- **Harmony** is a good default. It is fast, operates on the PCA
  embedding directly, and handles multiple batch variables
  simultaneously. It tends to work well when batch effects are moderate
  and cell type composition is broadly similar across batches.
  [`harmony_sc()`](https://gregorlueg.github.io/bixverse/reference/harmony_sc.md)
  runs version 2, which cares less about every cell type showing up in
  every batch and runs 4 k-means iterations per round instead of 20.
  [`harmony_v1_sc()`](https://gregorlueg.github.io/bixverse/reference/harmony_v1_sc.md)
  is only there to reproduce old results.

- **fastMNN** is worth considering when cell type composition differs
  substantially across batches. The mutual nearest neighbour approach is
  less sensitive to this because it only aligns cells that have a
  plausible biological match in another batch.

- **Seurat CCA and rPCA** are the anchor-based options. Like fastMNN
  they only align cells with a plausible match in another batch, so they
  hold up when cell type composition differs. They are the heaviest of
  the five though, and on data that already mixes well they can make the
  batch metrics slightly worse. rPCA is the cheaper and gentler of the
  two and a sensible default; go to CCA when the batches genuinely don’t
  share their cell types.

- **BBKNN** leaves the embedding alone and hands a batch-balanced graph
  straight to clustering. It’s the lightest touch of the lot, and that
  cuts both ways: when the batches sit far apart, as here, the
  cross-batch edges carry next to no weight and the clusters stay split
  by batch. Judge it on clustering, not on kBET or iLISI.

In practice, running multiple methods and looking at quantitative
metrics is the most reliable approach. No single metric captures
everything: kBET and PCR only see batch, cLISI and cell type ASW only
see biology, and the kNN metrics can’t score a method that builds the
kNN to look balanced. UMAPs and tSNE assessments should be taken with a
big pinch of salt, see [Chari, et
al.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011288)
Admittedly, they do look pretty however.

## Clean up

``` r

unlink(tempdir_batch_cor, recursive = TRUE, force = TRUE)
```
