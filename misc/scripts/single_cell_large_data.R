# large data set tests

devtools::load_all()

rextendr::clean()
rextendr::document()

tictoc::tic()
path_h5 <- path.expand(
  "~/Downloads/plate1_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad"
)

meta <- get_h5ad_dimensions(path_h5)

sc_object <- SingleCells(dir_data = path.expand("~/Desktop/bixverse_big_data/"))

sc_object <- stream_h5ad(object = sc_object, h5_path = path_h5)

# sc_object <- load_existing(object = sc_object)

var <- get_sc_var(object = sc_object)

gs_of_interest <- list(
  MT = var[grepl("^MT-", gene_id), gene_id],
  Ribo = var[grepl("^RPS|^RPL", gene_id), gene_id]
)

sc_object <- gene_set_proportions_sc(
  object = sc_object,
  gene_set_list = gs_of_interest,
  streaming = TRUE, # set streaming to TRUE !
  .verbose = TRUE
)

sc_object[[c(1:5L)]]

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

qc <- run_cell_qc(metrics, directions, threshold = 3)

qc

sc_object[["outlier"]] <- qc$combined

cells_to_keep <- qc_df[!qc$combined, cell_id]

sc_object <- set_cells_to_keep(sc_object, cells_to_keep)

sc_object <- find_hvg_sc(
  object = sc_object,
  streaming = TRUE # set streaming to TRUE!
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 30L,
  sparse_svd = TRUE # will automatically be set with ≥500k cells
)

# GPU-accelerated kNN search via a CAGRA style approach
# This will also add the shared nearest neighbour graph

library(bixverse.gpu)

# leverage a GPU-accelerated kNN search...
# otherwise, nndescent and HNSW are decent on the CPU
sc_object <- find_neighbours_cagra_sc(
  object = sc_object,
  cagra_params = params_sc_cagra(k_query = 15L, ann_dist = "cosine"),
  extract_knn = FALSE,
  .verbose = TRUE
)

# this will take time, but does actually run...
# but to be clear, this is running on 65m + edges
sc_object <- umap_sc(object = sc_object)

meta_data <- sc_object[[c("cell_id", "cell_name", "drug")]]

cells_grp_a <- meta_data[drug == "AT7519", cell_id]
cells_grp_b <- meta_data[drug == "palbociclib", cell_id]

dge_results <- find_markers_sc(
  object = sc_object,
  cells_1 = cells_grp_a,
  cells_2 = cells_grp_b
)

tictoc::toc()

setorder(dge_results, fdr)

head(dge_results, 25L)
