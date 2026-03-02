# large data set tests

devtools::load_all()

rextendr::clean()
rextendr::document()

path_h5 <- path.expand(
  "~/Downloads/plate1_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad"
)

sce <- single_cell_exp(dir_data = path.expand("~/Desktop/bixverse_big_data/"))

# sce <- stream_h5ad(object = sce, h5_path = path_h5)

sce <- load_existing(object = sce)

dim(sce)

pryr::object_size(sce)

tictoc::tic()
sce[1:100, , return_format = "cell"][1:5, 1:5]
tictoc::toc()

tictoc::tic()
sce[1:1000, , return_format = "cell", assay = "norm"][1:5, 1:5]
tictoc::toc()

tictoc::tic()
sce[, 1:10, return_format = "gene"][1:5, 1:5]
tictoc::toc()

sce <- find_hvg_sc(sce, streaming = TRUE)

sce <- calculate_pca_sc(sce, no_pcs = 30L)

pryr::object_size(sce)

sce <- find_neighbours_sc(
  object = sce,
  no_embd_to_use = 16L,
  neighbours_params = params_sc_neighbours(knn = list(knn_method = "hnsw"))
)

pryr::object_size(sce)

sce[1:100, , return_format = "cell"]

########
# HNSW #
########

# hnsw with 32 dim

# Building HNSW index with 2_857_393 nodes, M = 16
# Max layer: 6, Entry point: 299562
#   Layer 0: 2_857_393 nodes
#   Layer 1: 178_644 nodes
#   Layer 2: 11_036 nodes
#   Layer 3: 666 nodes
#   Layer 4: 33 nodes
#   Layer 5: 5 nodes
#   Layer 6: 1 nodes
# Building layer 6 with 1 nodes
#   Layer 6 built in 1.48ms
# Building layer 5 with 5 nodes
#   Layer 5 built in 1.72ms
# Building layer 4 with 33 nodes
#   Layer 4 built in 2.19ms
# Building layer 3 with 666 nodes
#   Layer 3 built in 12.58ms
# Building layer 2 with 11_036 nodes
#   Layer 2 built in 155.21ms
# Building layer 1 with 178_644 nodes
#   Layer 1 built in 4.31s
# Building layer 0 with 2_857_393 nodes
#   Layer 0 built in 106.59s
# Total HNSW build time: 111.38s
# Generated HNSW index: 111.42s
#   Processed 100_000 / 2_857_393 samples.
#   Processed 200_000 / 2_857_393 samples.
#   Processed 300_000 / 2_857_393 samples.
#   Processed 400_000 / 2_857_393 samples.
#   Processed 500_000 / 2_857_393 samples.
#   Processed 600_000 / 2_857_393 samples.
#   Processed 700_000 / 2_857_393 samples.
#   Processed 800_000 / 2_857_393 samples.
#   Processed 900_000 / 2_857_393 samples.
#   Processed 1_000_000 / 2_857_393 samples.
#   Processed 1_100_000 / 2_857_393 samples.
#   Processed 1_200_000 / 2_857_393 samples.
#   Processed 1_300_000 / 2_857_393 samples.
#   Processed 1_400_000 / 2_857_393 samples.
#   Processed 1_500_000 / 2_857_393 samples.
#   Processed 1_600_000 / 2_857_393 samples.
#   Processed 1_700_000 / 2_857_393 samples.
#   Processed 1_800_000 / 2_857_393 samples.
#   Processed 1_900_000 / 2_857_393 samples.
#   Processed 2_000_000 / 2_857_393 samples.
#   Processed 2_100_000 / 2_857_393 samples.
#   Processed 2_200_000 / 2_857_393 samples.
#   Processed 2_300_000 / 2_857_393 samples.
#   Processed 2_400_000 / 2_857_393 samples.
#   Processed 2_500_000 / 2_857_393 samples.
#   Processed 2_600_000 / 2_857_393 samples.
#   Processed 2_700_000 / 2_857_393 samples.
#   Processed 2_800_000 / 2_857_393 samples.
# Identified approximate nearest neighbours via HNSW: 59.52s
# Recall of approximate nearest neighbours search in random subset: 1.00
# KNN generation done : 185.52s
# Generating sNN graph (full: FALSE).
# Transformed kNN into an sNN graph: 13.64s
# Transforming sNN data to igraph.

#########
# Annoy #
#########

# annoy with 32 dims

# Generated Annoy index: 15.24s
#   Processed 100_000 / 2_857_393 samples.
#   Processed 200_000 / 2_857_393 samples.
#   Processed 300_000 / 2_857_393 samples.
#   Processed 400_000 / 2_857_393 samples.
#   Processed 500_000 / 2_857_393 samples.
#   Processed 600_000 / 2_857_393 samples.
#   Processed 700_000 / 2_857_393 samples.
#   Processed 800_000 / 2_857_393 samples.
#   Processed 900_000 / 2_857_393 samples.
#   Processed 1_000_000 / 2_857_393 samples.
#   Processed 1_100_000 / 2_857_393 samples.
#   Processed 1_200_000 / 2_857_393 samples.
#   Processed 1_300_000 / 2_857_393 samples.
#   Processed 1_400_000 / 2_857_393 samples.
#   Processed 1_500_000 / 2_857_393 samples.
#   Processed 1_600_000 / 2_857_393 samples.
#   Processed 1_700_000 / 2_857_393 samples.
#   Processed 1_800_000 / 2_857_393 samples.
#   Processed 1_900_000 / 2_857_393 samples.
#   Processed 2_000_000 / 2_857_393 samples.
#   Processed 2_100_000 / 2_857_393 samples.
#   Processed 2_200_000 / 2_857_393 samples.
#   Processed 2_300_000 / 2_857_393 samples.
#   Processed 2_400_000 / 2_857_393 samples.
#   Processed 2_500_000 / 2_857_393 samples.
#   Processed 2_600_000 / 2_857_393 samples.
#   Processed 2_700_000 / 2_857_393 samples.
#   Processed 2_800_000 / 2_857_393 samples.
# Identified approximate nearest neighbours via Annoy: 371.87s
# Recall of approximate nearest neighbours search in random subset: 0.99
# KNN generation done : 402.73s
# Generating sNN graph (full: FALSE).
# Transformed kNN into an sNN graph: 14.19s
# Transforming sNN data to igraph.

#############
# NNDescent #
#############

# Built Annoy index: 9.18s
# Running NN-Descent: 2_857_393 samples, max_candidates=30
# Queried Annoy index: 53.26s
# Using chunk size 145_635 (20 chunks) for memory-efficient updates
#  Preparing candidates for iter 1
#  Processing updates for iter 1 (20 chunks)
#   Iter 1: 41_134_731 edge updates (rate=0.4799)
#  Preparing candidates for iter 2
#  Processing updates for iter 2 (20 chunks)
#   Iter 2: 8_975_855 edge updates (rate=0.1047)
#  Preparing candidates for iter 3
#  Processing updates for iter 3 (20 chunks)
#   Iter 3: 908_498 edge updates (rate=0.0106)
#  Preparing candidates for iter 4
#  Processing updates for iter 4 (20 chunks)
#   Iter 4: 63_648 edge updates (rate=0.0007)
#   Converged after 4 iterations
# Total time: 204.08s
# Generated NNDescent index: 214.81s
#   Processed 100_000 / 2_857_393 samples
#   Processed 200_000 / 2_857_393 samples
#   Processed 300_000 / 2_857_393 samples
#   Processed 400_000 / 2_857_393 samples
#   Processed 500_000 / 2_857_393 samples
#   Processed 600_000 / 2_857_393 samples
#   Processed 700_000 / 2_857_393 samples
#   Processed 800_000 / 2_857_393 samples
#   Processed 900_000 / 2_857_393 samples
#   Processed 1_000_000 / 2_857_393 samples
#   Processed 1_100_000 / 2_857_393 samples
#   Processed 1_200_000 / 2_857_393 samples
#   Processed 1_300_000 / 2_857_393 samples
#   Processed 1_400_000 / 2_857_393 samples
#   Processed 1_500_000 / 2_857_393 samples
#   Processed 1_600_000 / 2_857_393 samples
#   Processed 1_700_000 / 2_857_393 samples
#   Processed 1_800_000 / 2_857_393 samples
#   Processed 1_900_000 / 2_857_393 samples
#   Processed 2_000_000 / 2_857_393 samples
#   Processed 2_100_000 / 2_857_393 samples
#   Processed 2_200_000 / 2_857_393 samples
#   Processed 2_300_000 / 2_857_393 samples
#   Processed 2_400_000 / 2_857_393 samples
#   Processed 2_500_000 / 2_857_393 samples
#   Processed 2_600_000 / 2_857_393 samples
#   Processed 2_700_000 / 2_857_393 samples
#   Processed 2_800_000 / 2_857_393 samples
# Identified approximate nearest neighbours via NNDescent: 17.70s.
# Recall of approximate nearest neighbours search in random subset: 1.00
# KNN generation done : 246.43s
# Generating sNN graph (full: FALSE).
# Transformed kNN into an sNN graph: 17.87s
# Transforming sNN data to igraph.
