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
  neighbours across batches – pairs of cells from different batches that
  are each other’s closest match. Correction vectors are computed from
  these pairs and applied to align the embeddings. It is particularly
  effective when cell type composition differs across batches.

- **Harmony** ([Korsunsky, et al.,
  2019](https://doi.org/10.1038/s41592-019-0619-0)) operates on PCA
  embeddings using iterative soft clustering with diversity penalties.
  It assigns cells to clusters, estimates batch effects per cluster via
  ridge regression, and corrects the embedding. It tends to be fast and
  works well across a range of scenarios. Additionally, a version 2 has
  been released, see [Patikas, et al.,
  2026](https://www.biorxiv.org/content/10.64898/2026.03.16.711825v1)
  with improvements that increase scalability to larger data sets and
  modifications to the objective functions to not overcorrect across
  batches when different cell types are there.

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
#>  Using light streaming for the CSR to CSC conversion.
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
#> Using dense SVD solving on scaled data on 2000 HVG.

sc_object <- find_neighbours_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours()
)
#> 
#> Generating sNN graph (full: TRUE).
#> Transforming sNN data to igraph.
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
#>       874        79      1016       331        12      3529
```

## Comparing the different methods

Let’s check out the different methods and how they behave.

### Uncorrected data

Before applying any correction, we can visualise the batch effect and
compute baseline metrics.

``` r

sc_object <- tsne_sc(object = sc_object)
#> Running t-SNE.
```

tSNE plot (tSNE’s usually look prettier… It’s anyway not a good way to
assess batch effect correction, so at least let’s make it pleasant on
the eye).

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
#>   Cells: 5841 | Batches: 2 | Threshold: 0.050
#>   Rejection rate:      0.9885 (5774 / 5841)
#>   Mean Chi-Square:     12.8259 (expected under H0: 1)
#>   Median Chi-Square:   8.7578
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
#> 1:       pca  0.01147064 0.8909914     0             NA     1     0.6758695
#>    graph_connectivity      method
#>                 <num>      <char>
#> 1:          0.9894515 Uncorrected
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
#> Generating sNN graph (full: TRUE).
#> Transforming sNN data to igraph.
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
#> 1:       mnn   0.6449238 0.9378912 0.3005779      0.8858816     1     0.6544492
#>    graph_connectivity  method
#>                 <num>  <char>
#> 1:          0.9972923 fastMNN
```

We can see clear improvements on the batch side: kBET acceptance, batch
ASW and iLISI all go up.

``` r

sc_object <- tsne_sc(object = sc_object, slot_name = "tsne_mnn", use_knn = TRUE)
#> Running t-SNE.
#> Using provided kNN graph.
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

Harmony operates on the PCA embedding directly. The number of clusters
is auto-determined from the dataset size (capped at 100) when left as
`NULL`.

The two versions of Harmony are supported. Let’s start with the
original:

#### Version 1

``` r

sc_object <- harmony_sc(
  object = sc_object,
  batch_column = "exp_id",
  harmony_params = params_sc_harmony()
)
#>  Auto-determined number of Harmony clusters: 100

sc_object <- find_neighbours_sc(
  object = sc_object,
  embd_to_use = "harmony",
  neighbours_params = params_sc_neighbours()
)
#> 
#> Generating sNN graph (full: TRUE).
#> Transforming sNN data to igraph.
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
    .[, method := "Harmony"]
)

metrics_dt[method == "Harmony"]
#>    embedding kbet_accept batch_asw     ilisi pcr_comparison clisi cell_type_asw
#>       <char>       <num>     <num>     <num>          <num> <num>         <num>
#> 1:   harmony   0.8635508 0.9270982 0.6423357      0.9235815     1      0.647598
#>    graph_connectivity  method
#>                 <num>  <char>
#> 1:           0.989357 Harmony
```

Also here, we observe improvements across the board.

``` r

sc_object <- tsne_sc(
  object = sc_object,
  slot_name = "tsne_harm",
  use_knn = TRUE
)
#> Running t-SNE.
#> Using provided kNN graph.
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
    title = "Harmony batch correction",
    colour = "Batch:"
  )
```

![](single_cell_batch_corrections_files/figure-html/harmony%20plot-1.png)

#### Version 2

Compared to v1 of Harmony, several modifications/improvements were
implemented by [Patikas, et al.,
2026](https://www.biorxiv.org/content/10.64898/2026.03.16.711825v1),
namely:

- Stabilised diversity penalty
- Batch pruning in ridge regression
- Arrowhead matrix inversion for single-covariate case (making it
  faster!)
- Dynamic lambda estimation and theta scaling by batch size.

Overall, this makes v2 more amenable and faster on big data sets than v1
(if you wish to use Harmony, maybe use the version 2 over the v2?).

``` r

sc_object <- harmony_v2_sc(
  object = sc_object,
  batch_column = "exp_id",
  harmony_params = params_sc_harmony_v2()
)
#>  Auto-determined number of Harmony clusters: 100

sc_object <- find_neighbours_sc(
  object = sc_object,
  embd_to_use = "harmony_v2",
  neighbours_params = params_sc_neighbours()
)
#> 
#> Generating sNN graph (full: TRUE).
#> Transforming sNN data to igraph.
```

Let’s calculate the batch correction-related scores

``` r

metrics_dt <- rbind(
  metrics_dt,
  calculate_integration_metrics_sc(
    sc_object,
    batch_column = "exp_id",
    cell_type_column = "cell_type",
    embd_to_use = "harmony_v2",
    .verbose = FALSE
  ) %>%
    .[, method := "Harmony v2"]
)

metrics_dt[method == "Harmony v2"]
#>     embedding kbet_accept batch_asw     ilisi pcr_comparison clisi
#>        <char>       <num>     <num>     <num>          <num> <num>
#> 1: harmony_v2   0.8604691 0.9189611 0.6423357      0.8947806     1
#>    cell_type_asw graph_connectivity     method
#>            <num>              <num>     <char>
#> 1:     0.6957046          0.9879408 Harmony v2
```

``` r

sc_object <- tsne_sc(
  object = sc_object,
  slot_name = "tsne_harm_v2",
  use_knn = TRUE
)
#> Running t-SNE.
#> Using provided kNN graph.
```

``` r

embedding_plot_sc(
  sc_object,
  embedding = "tsne_harm_v2",
  colour_by = "exp_id",
  label_by = "exp_id",
  discrete = TRUE
) +
  labs(
    title = "Harmony (version 2) batch correction",
    colour = "Batch:"
  )
```

![](single_cell_batch_corrections_files/figure-html/harmony%20(v2)%20plot-1.png)

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
#> Generating sNN graph (full: TRUE).
#> Transforming sNN data to igraph.
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
#> 1:       cca   0.6856703 0.8886744 0.4705881      0.8693617     1     0.6997823
#>    graph_connectivity     method
#>                 <num>     <char>
#> 1:          0.9894515 Seurat CCA
```

``` r

sc_object <- tsne_sc(object = sc_object, slot_name = "tsne_cca", use_knn = TRUE)
#> Running t-SNE.
#> Using provided kNN graph.
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
#> Generating sNN graph (full: TRUE).
#> Transforming sNN data to igraph.
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
#> 1:      rpca   0.7752097 0.9051192 0.6423357       0.810686     1      0.691866
#>    graph_connectivity      method
#>                 <num>      <char>
#> 1:          0.9894043 Seurat rPCA
```

``` r

sc_object <- tsne_sc(
  object = sc_object,
  slot_name = "tsne_rpca",
  use_knn = TRUE
)
#> Running t-SNE.
#> Using provided kNN graph.
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

BBKNN constructs a batch-balanced kNN graph rather than correcting an
embedding. The key parameter is `neighbours_within_batch`, which
controls how many neighbours are found per batch. Setting
`no_neighbours_to_keep` below the total generated neighbours (here 10
from 2 batches = 20 total) ensures the distance-based filtering
introduces genuine variation, making downstream metrics informative.

Note that kBET and the embedding metrics are not well suited for
evaluating BBKNN: kBET compares neighbourhood proportions against global
proportions, which BBKNN manipulates by design, and ASW and PCR need an
embedding. BBKNN only hands back a graph, so we pass
`embd_to_use = NULL` and those columns come back `NA`. iLISI, cLISI and
graph connectivity on the stored kNN are the ones to read here.

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
#>    embedding kbet_accept batch_asw ilisi pcr_comparison clisi cell_type_asw
#>       <char>       <num>     <num> <num>          <num> <num>         <num>
#> 1:      <NA>           1        NA   0.8             NA     1            NA
#>    graph_connectivity method
#>                 <num> <char>
#> 1:          0.9755626  BBKNN
```

``` r

sc_object <- tsne_sc(object = sc_object, use_knn = TRUE)
#> Running t-SNE.
#> Using provided kNN graph.
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
#> 1: Uncorrected  0.01147064 0.8909914 0.0000000             NA     1
#> 2:     fastMNN  0.64492381 0.9378912 0.3005779      0.8858816     1
#> 3:     Harmony  0.86355076 0.9270982 0.6423357      0.9235815     1
#> 4:  Harmony v2  0.86046910 0.9189611 0.6423357      0.8947806     1
#> 5:  Seurat CCA  0.68567026 0.8886744 0.4705881      0.8693617     1
#> 6: Seurat rPCA  0.77520972 0.9051192 0.6423357      0.8106860     1
#> 7:       BBKNN  1.00000000        NA 0.8000000             NA     1
#>    cell_type_asw graph_connectivity
#>            <num>              <num>
#> 1:     0.6758695          0.9894515
#> 2:     0.6544492          0.9972923
#> 3:     0.6475980          0.9893570
#> 4:     0.6957046          0.9879408
#> 5:     0.6997823          0.9894515
#> 6:     0.6918660          0.9894043
#> 7:            NA          0.9755626
```

``` r

bio_metrics <- c("clisi", "cell_type_asw", "graph_connectivity")

metrics_long <- melt(
  metrics_dt[, -"embedding"],
  id.vars = "method",
  variable.name = "metric",
  value.name = "score"
) %>%
  .[, `:=`(
    method = factor(method, levels = rev(unique(metrics_dt$method))),
    type = fifelse(metric %in% bio_metrics, "Bio conservation", "Batch mixing")
  )]

ggplot(metrics_long, aes(x = metric, y = method, fill = score)) +
  geom_tile(colour = "white") +
  geom_text(
    aes(label = fifelse(is.na(score), "", sprintf("%.2f", score))),
    size = 3
  ) +
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
  labs(x = NULL, y = NULL, fill = "Score", title = "Integration metrics")
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
two PBMC runs, kBET, iLISI and PCR do the talking. And BBKNN’s perfect
kBET is the artefact mentioned above, not a win.

There is no universally best batch correction method. Some practical
guidance:

- **Harmony** is a good default. It is fast, operates on the PCA
  embedding directly, and handles multiple batch variables
  simultaneously. It tends to work well when batch effects are moderate
  and cell type composition is broadly similar across batches. Version 2
  is quite an improvement and necessitates less that the same cells are
  available across the batches and is substantially faster when you
  provide only one batch effect to regress out.

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

- **BBKNN** is useful when you want to avoid modifying the embedding at
  all. The batch-balanced graph can be passed directly to graph-based
  clustering algorithms. It is the lightest-touch approach but provides
  less control over the degree of correction.

In practice, running multiple methods and looking at quantitative
metrics is the most reliable approach. No single metric captures
everything: kBET and PCR only see batch, cLISI and cell type ASW only
see biology, and the ones built on the kNN graph are the only fair
comparison for graph-based methods like BBKNN. UMAPs and tSNE
assessments should be taken with a big pinch of salt, see [Chari, et
al.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011288)
Admittedly, they do look pretty however.

## Clean up

``` r

unlink(tempdir_batch_cor, recursive = TRUE, force = TRUE)
```
