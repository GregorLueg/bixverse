# large data set tests

devtools::load_all()

rextendr::clean()
rextendr::document()

path_h5 <- path.expand(
  "~/Downloads/plate1_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad"
)

sce <- single_cell_exp(dir_data = path.expand("~/Desktop/bixverse_big_data/"))

sce <- load_existing(object = sce)

# sce <- stream_h5ad(object = sce, h5_path = path_h5)

sce <- find_hvg_sc(sce, streaming = TRUE)

sce <- calculate_pca_sc(sce, no_pcs = 32L)

sce <- find_neighbours_sc(
  object = sce,
  no_embd_to_use = 15L,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "ivf", n_centroids = 1000L, n_probes = 15L)
  )
)

########
# HNSW #
########

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
#   Layer 6 built in 2.13ms
# Building layer 5 with 5 nodes
#   Layer 5 built in 3.73ms
# Building layer 4 with 33 nodes
#   Layer 4 built in 4.66ms
# Building layer 3 with 666 nodes
#   Layer 3 built in 10.28ms
# Building layer 2 with 11_036 nodes
#   Layer 2 built in 143.95ms
# Building layer 1 with 178_644 nodes
#   Layer 1 built in 4.58s
# Building layer 0 with 2_857_393 nodes
#   Layer 0 built in 92.02s
# Total HNSW build time: 97.06s
# Generated HNSW index: 97.08s
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
# Identified approximate nearest neighbours via HNSW: 52.29s
# Recall of approximate nearest neighbours search in random subset: 1.00
# KNN generation done : 165.77s
# Generating sNN graph (full: FALSE).
# Transformed kNN into an sNN graph: 15.17s
# Transforming sNN data to igraph.

#########
# Annoy #
#########

# Generated Annoy index: 11.66s
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
# Identified approximate nearest neighbours via Annoy: 121.47s
# Recall of approximate nearest neighbours search in random subset: 1.00
# KNN generation done : 149.60s
# Generating sNN graph (full: FALSE).
# Transformed kNN into an sNN graph: 15.56s
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

#######
# IVF #
#######

# with n_centroids = 1000L
# and n_probe = 15L

#  Sampling 250k vectors for training
#   Initialising centroids via fast random selection
#   Running parallel Lloyd's iterations
#     Iteration 1 complete
#     Iteration 2 complete
#     Iteration 3 complete
#     Iteration 4 complete
#     Iteration 5 complete
#     Iteration 6 complete
#     Iteration 7 complete
#     Iteration 8 complete
#     Iteration 9 complete
#     Iteration 10 complete
#     Iteration 11 complete
#     Iteration 12 complete
#     Iteration 13 complete
#     Iteration 14 complete
#     Iteration 15 complete
#     Iteration 16 complete
#     Iteration 17 complete
#     Iteration 18 complete
#     Iteration 19 complete
#     Iteration 20 complete
#     Iteration 21 complete
#     Iteration 22 complete
#     Iteration 23 complete
#     Iteration 24 complete
#     Iteration 25 complete
#     Iteration 26 complete
#     Iteration 27 complete
#     Iteration 28 complete
#     Iteration 29 complete
#     Iteration 30 complete
#     Iteration 31 complete
#     Iteration 32 complete
#     Iteration 33 complete
#     Iteration 34 complete
#     Iteration 35 complete
#     Iteration 36 complete
#     Iteration 37 complete
#     Iteration 38 complete
#     Iteration 39 complete
#     Iteration 40 complete
#     Iteration 41 complete
#     Iteration 42 complete
#     Iteration 43 complete
#     Iteration 44 complete
#     Iteration 45 complete
#     Iteration 46 complete
#     Iteration 47 complete
#     Iteration 48 complete
#     Iteration 49 complete
#     Iteration 50 complete
# Generated IVF index: 5.74s
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
# Identified approximate nearest neighbours via IVF index: 768.89s.
# Recall of approximate nearest neighbours search in random subset: 1.00
# KNN generation done : 801.62s
# Generating sNN graph (full: FALSE).
# Transformed kNN into an sNN graph: 16.17s
# Transforming sNN data to igraph.
