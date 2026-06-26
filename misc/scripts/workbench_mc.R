library(bixverse)
library(ggplot2)
library(data.table)
library(bixverse.plots)
library(magrittr)
pbmc3k_path <- download_pbmc3k()

tempdir_pbmc <- tempdir()

sc_object <- SingleCells(dir_data = tempdir_pbmc)

mtx_io_params <- get_cell_ranger_params(pbmc3k_path)


sc_object <- load_mtx(
  object = sc_object,
  sc_mtx_io_param = mtx_io_params,
  mtx_streaming = FALSE,
  .verbose = TRUE
)
sc_object[["sample"]] <- sample(
  x = c("sample1", "sample2"),
  size = nrow(sc_object),
  replace = T
)

grp_object <- SingleCellsSubset(
  sc_object = sc_object,
  grouping_column = "sample",
  group = "sample1"
)
head(grp_object)
print(grp_object)
get_sc_var(grp_object)
get_sc_counts(grp_object)


object = grp_object
group_by = "sample"
sc_meta_cell_params = params_sc_bt_metacells()
regenerate_knn = FALSE
embd_to_use = "pca"
no_embd_to_use = NULL
cells_to_use = NULL
target_size = 1e5
seed = 42L
.verbose = TRUE


##
setnames_sc(
  object = sc_object,
  table = "var",
  old = "column1",
  new = "gene_symbol"
)

var <- get_sc_var(sc_object)

ensembl_to_symbol <- setNames(var$gene_symbol, var$gene_id)
symbol_to_ensembl <- setNames(var$gene_id, var$gene_symbol)

# Define gene sets whose proportional expression signals low quality
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

# Collect metrics for outlier detection
qc_df <- sc_object[[c("cell_id", "sample", "lib_size", "nnz", "MT")]]

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
  threshold = 3,
  groups = qc_df$sample
)

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = 2000L,
  .verbose = TRUE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 30L,
  sparse_svd = TRUE
)

# the data is so tiny that exhaustive kNN search is faster than building
# an approximate nearest neighbour index
sc_object <- find_neighbours_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "exhaustive")
  )
)
sc_object <- find_clusters_sc(sc_object, res = 1, name = "leiden_clusters")
sc_object <- umap_sc(sc_object)

cell_markers <- c(
  CD3D = "T cells",
  CD3E = "T cells",
  CD3G = "T cells",
  IL7R = "CD4+ T",
  CD4 = "CD4+ T",
  CD8A = "CD8+ T",
  CD8B = "CD8+ T",
  MS4A1 = "B cells",
  CD79A = "B cells",
  CD19 = "B cells",
  CD14 = "CD14+ Mono",
  LYZ = "CD14+ Mono",
  S100A8 = "CD14+ Mono",
  FCGR3A = "CD16+ Mono",
  CDKN1C = "CD16+ Mono",
  GNLY = "NK",
  NKG7 = "NK",
  NCAM1 = "NK",
  FCER1A = "mDC",
  CD1C = "mDC",
  LILRA4 = "pDC",
  CLEC4C = "pDC",
  PPBP = "Platelet",
  PF4 = "Platelet"
)
cell_markers_dt <- stack(cell_markers) %>%
  as.data.table() %>%
  setnames(., c("values", "ind"), c("cell_type", "gene_symbol")) %>%
  .[, gene_symbol := as.character(gene_symbol)] %>%
  .[, gene_id := var$gene_id[match(gene_symbol, var$gene_symbol)]] %>%
  .[!is.na(gene_id), ]

## Prepare the list of markers
marker_list <- prepare_cell_markers(sc_object, cell_markers_dt)
## Calculate the score for each cell
sctype_scores <- calc_sc_type_scores(
  object = sc_object,
  cell_marker_list = marker_list
)
## Annotate the clusters
cell_type_anno <- score_clusters(
  sctype_scores,
  sc_object[[]][["leiden_clusters"]]
)

obs <- sc_object[[]][, .(cell_idx, leiden_clusters)]
add_sctype <- setNames(
  cell_type_anno$cell_type[match(
    obs$leiden_clusters,
    cell_type_anno$cluster_id
  )],
  nm = obs$cell_idx
)
sc_object[["sc_type"]] <- add_sctype


##
mc <- generate_supercells_sc(
  object = sc_object,
  .verbose = TRUE
)

sc_supercell_params = params_sc_supercell()
regenerate_knn = FALSE
embd_to_use = "pca"
no_embd_to_use = NULL
group = NULL
cells_to_use = NULL
target_size = 1e5
seed = 42L
.verbose = TRUE
