# large data set tests

path_h5 <- path.expand(
  "~/Downloads/plate1_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad"
)

sce <- single_cell_exp(dir_data = tempdir())

sce <- stream_h5ad(object = sce, h5_path = path_h5)

sce <- find_hvg_sc(sce, streaming = TRUE)

sce <- calculate_pca_sc(sce, no_pcs = 30L)

sce <- find_neighbours_sc(
  object = sce,
  no_embd_to_use = 15L,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "annoy")
  )
)

## Annoy speed

# Generated Annoy index: 31.88s
#  Processed 100_000 / 2_857_393 samples.
#  Processed 200_000 / 2_857_393 samples.
#  Processed 300_000 / 2_857_393 samples.
#  Processed 400_000 / 2_857_393 samples.
#  Processed 500_000 / 2_857_393 samples.
#  Processed 600_000 / 2_857_393 samples.
#  Processed 700_000 / 2_857_393 samples.
#  Processed 800_000 / 2_857_393 samples.
#  Processed 900_000 / 2_857_393 samples.
#  Processed 1_000_000 / 2_857_393 samples.
#  Processed 1_100_000 / 2_857_393 samples.
#  Processed 1_200_000 / 2_857_393 samples.
#  Processed 1_300_000 / 2_857_393 samples.
#  Processed 1_400_000 / 2_857_393 samples.
#  Processed 1_500_000 / 2_857_393 samples.
#  Processed 1_600_000 / 2_857_393 samples.
#  Processed 1_700_000 / 2_857_393 samples.
#  Processed 1_800_000 / 2_857_393 samples.
#  Processed 1_900_000 / 2_857_393 samples.
#  Processed 2_000_000 / 2_857_393 samples.
#  Processed 2_100_000 / 2_857_393 samples.
#  Processed 2_200_000 / 2_857_393 samples.
#  Processed 2_300_000 / 2_857_393 samples.
#  Processed 2_400_000 / 2_857_393 samples.
#  Processed 2_500_000 / 2_857_393 samples.
#  Processed 2_600_000 / 2_857_393 samples.
#  Processed 2_700_000 / 2_857_393 samples.
#  Processed 2_800_000 / 2_857_393 samples.
# Identified approximate nearest neighbours via Annoy: 129.39s.
# KNN generation done : 167.21s
# Generating sNN graph (full: FALSE).
# Transformed kNN into an sNN graph: 16.88s
# Transforming sNN data to igraph.

## HNSW

# Building HNSW index with 2_857_393 nodes, M = 32
# Max layer: 4, Entry point: 156061
# Building layer 4 with 5 nodes
#   Layer 4 built in 6.10ms
# Building layer 3 with 71 nodes
#   Layer 3 built in 3.57ms
# Building layer 2 with 2_747 nodes
#   Layer 2 built in 55.44ms
# Building layer 1 with 89_236 nodes
#   Layer 1 built in 3.18s
# Building layer 0 with 2_857_393 nodes
#   Layer 0 built in 161.00s
# Total HNSW build time: 164.67s
# Generated HNSW index: 164.67s
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
# Identified approximate nearest neighbours via HNSW: 86.13s.
# KNN generation done : 250.81s
# Generating sNN graph (full: FALSE).
# Transformed kNN into an sNN graph: 16.82s
# Transforming sNN data to igraph.
