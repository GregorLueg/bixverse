# atac against seurat -----------------------------------------------------------

# port of bixverse-rs python/tests/test_atac_seurat.py: LSI projection, rLSI
# integration and label transfer against Seurat 6.0.0. Where an SVD sign would
# differ, Seurat's own LSI goes in.

source("helper_atac.R", local = TRUE)
if (!ATAC_OK) {
  exit_file("ATAC fixtures not found; set BIXVERSE_ATAC_FIXTURES")
}

batch <- as.integer(atac_lines("seurat/batch.txt"))
lsi_features <- as.integer(atac_lines("seurat/features.txt"))
exact <- list(knn_method = "exhaustive")

tfidf_experiment <- function(barcodes) {
  run_tfidf_sc(atac_mm_experiment(barcodes))
}

## project lsi -----------------------------------------------------------------

ref <- tfidf_experiment(ATAC_BARCODES[batch == 0L])
query <- tfidf_experiment(ATAC_BARCODES[batch == 1L])
ref <- calculate_lsi_sc(
  ref,
  features = lsi_features + 1L,
  lsi_params = params_atac_lsi(
    n_comp = 30L,
    drop_first = FALSE,
    solver = "exact"
  ),
  .verbose = FALSE
)
info <- ref@atac_info
info$lsi_loadings <- atac_dense("seurat/batch0_lsi_loadings.tsv.gz")
info$lsi_singular_values <- atac_tsv("seurat/batch0_lsi_d.tsv")$d
ref@atac_info <- info

query <- project_lsi_sc(ref, query)
got <- get_embedding(query, "lsi_projected", modality = "atac")
expect_equal(
  unname(got),
  unname(atac_dense("seurat/project_batch1_on_batch0.tsv.gz")),
  tolerance = 1e-6,
  info = "projection matches Seurat"
)

ref5 <- calculate_lsi_sc(
  tfidf_experiment(ATAC_BARCODES[1:250]),
  lsi_params = params_atac_lsi(n_comp = 5L, solver = "exact"),
  .verbose = FALSE
)
tiles <- add_atac_fragments_sc(
  atac_mm_experiment(ATAC_BARCODES[251:500], peaks = FALSE),
  ATAC_FRAGS,
  tile_width = 5000L,
  chrom_sizes = c(chr22 = 50818468L),
  mode = "fragment",
  .verbose = FALSE
)
tiles <- run_tfidf_sc(tiles)
expect_error(
  project_lsi_sc(ref5, tiles),
  pattern = "not in the query",
  info = "the query needs the reference features"
)

## rlsi ------------------------------------------------------------------------

exp <- tfidf_experiment(ATAC_BARCODES)
get_sc_duckdb(exp)$add_data_obs(data.table::data.table(seurat_batch = batch))
exp <- calculate_lsi_sc(
  exp,
  features = lsi_features + 1L,
  lsi_params = params_atac_lsi(
    n_comp = 30L,
    drop_first = FALSE,
    solver = "exact"
  ),
  .verbose = FALSE
)
# Seurat's combined LSI, so the comparison does not hinge on SVD signs
exp <- set_embedding(
  exp,
  embd = atac_dense("seurat/combined_lsi_embeddings.tsv.gz"),
  name = "lsi",
  modality = "atac"
)
exp <- run_rlsi_sc(
  exp,
  "seurat_batch",
  rlsi_params = params_atac_rlsi(
    k_weight = 20L,
    solver = "exact",
    knn = exact
  ),
  .verbose = FALSE
)
expect_equal(
  unname(get_embedding(exp, "rlsi", modality = "atac")),
  unname(atac_dense("seurat/rlsi_integrated.tsv.gz")),
  tolerance = 1e-5,
  info = "rLSI matches Seurat"
)

## label transfer --------------------------------------------------------------

# Ensembl IDs without version, as Cell Ranger keys the RNA
genes <- sub("\\..*$", "", atac_tsv("signac/gene_activity/genes.tsv")$gene_id)
rna <- atac_tsv("seurat/rna/matrix.tsv.gz")
counts <- Matrix::sparseMatrix(
  i = rna$cell_idx + 1L,
  j = rna$feature_idx + 1L,
  x = as.numeric(rna$value),
  dims = c(length(ATAC_BARCODES), length(genes)),
  dimnames = list(ATAC_BARCODES, genes)
)
counts <- as(counts, "RsparseMatrix")
ref_cells <- which(batch == 0L)
query_cells <- which(batch == 1L)
labels <- rep(NA_character_, length(ATAC_BARCODES))
labels[ref_cells] <- atac_lines("seurat/transfer_labels.txt")

dir_rna <- tempfile("atac_rna")
dir.create(dir_rna)
rna_ref <- load_r_data(
  SingleCells(dir_rna),
  counts = counts,
  obs = data.table::data.table(cell_id = ATAC_BARCODES, cell_type = labels),
  var = data.table::data.table(gene_id = genes),
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 0L,
    min_lib_size = 0L,
    min_cells = 0L,
    target_size = 1e4
  ),
  .verbose = FALSE
)
expect_equal(dim(rna_ref), c(length(ATAC_BARCODES), length(genes)))
rna_ref <- set_cells_to_keep(rna_ref, ATAC_BARCODES[ref_cells])

target <- as.numeric(atac_lines("seurat/transfer_ga_target.txt")[1])
ga <- gene_activity_sc(
  atac_mm_experiment(),
  ATAC_FRAGS,
  ATAC_GTF,
  target_size = target,
  .verbose = FALSE
)
ga <- set_cells_to_keep(ga, ATAC_BARCODES[query_cells])

features <- genes[as.integer(atac_lines("seurat/transfer_features.txt")) + 1L]
lsi <- atac_dense("seurat/batch1_lsi_embeddings.tsv.gz")
out <- transfer_labels_seurat_sc(
  rna_ref,
  ga,
  "cell_type",
  weight_embedding = lsi[, 2:30],
  features = features,
  transfer_params = params_sc_seurat_transfer(
    k_weight = as.integer(atac_lines("seurat/transfer_k_weight.txt")[1]),
    knn = exact
  ),
  .verbose = FALSE
)
want <- atac_tsv("seurat/transfer_predictions_exact_knn.tsv")
score_cols <- setdiff(names(want), names(want)[1])
lab_names <- sub("^prediction\\.score\\.", "", score_cols)
ours <- as.matrix(out[, sprintf("score_%s", lab_names), with = FALSE])
theirs <- as.matrix(want[, score_cols, with = FALSE])
diff <- apply(abs(ours - theirs), 1L, max)
# anchors swap between tied zero-count reference cells, which moves the
# scores of the query cells they weight; which ties flip is down to GEMM
# rounding, so a few cells move
expect_true(max(diff) < 1e-2, info = "scores match Seurat")
expect_true(sum(diff > 1e-5) <= 75L, info = "few cells move")
expect_equal(
  out$predicted_cell_type,
  lab_names[max.col(theirs, ties.method = "first")],
  info = "predictions match Seurat"
)
expect_equal(
  get_sc_obs(ga, cols = "predicted_cell_type", filtered = TRUE)[[1]],
  out$predicted_cell_type,
  info = "predictions written to the query obs"
)
