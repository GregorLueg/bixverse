# test scrublet ----------------------------------------------------------------

## parameters ------------------------------------------------------------------

n_doublets <- 200

test_wide_scrublet_params <- params_scrublet(
  no_pcs = 5L,
  expected_doublet_rate = 0.2,
  sim_doublet_ratio = 2.0,
  min_gene_var_pctl = 0.7,
  k = 3L
)

## helper functions ------------------------------------------------------------

#' Helper function to calculate metrics
metrics_helper <- function(cm) {
  TP <- cm[2, 2]
  FP <- cm[1, 2]
  FN <- cm[2, 1]
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1 <- 2 * (precision * recall) / (precision + recall)
  c(precision = precision, recall = recall, f1 = f1)
}

## initial data ----------------------------------------------------------------

syn_data <- generate_single_cell_test_data()

# Get cell indices by cell type
ct1_idx <- which(syn_data$obs$cell_grp == "cell_type_1")
ct2_idx <- which(syn_data$obs$cell_grp == "cell_type_2")
ct3_idx <- which(syn_data$obs$cell_grp == "cell_type_3")

# Create cross-cell-type doublets (more detectable)
n_doublets_12 <- ceiling(n_doublets / 3)
n_doublets_23 <- ceiling(n_doublets / 3)
n_doublets_13 <- n_doublets - n_doublets_12 - n_doublets_23

doublet_pairs <- rbind(
  cbind(
    sample(ct1_idx, n_doublets_12, replace = TRUE),
    sample(ct2_idx, n_doublets_12, replace = TRUE)
  ),
  cbind(
    sample(ct2_idx, n_doublets_23, replace = TRUE),
    sample(ct3_idx, n_doublets_23, replace = TRUE)
  ),
  cbind(
    sample(ct1_idx, n_doublets_13, replace = TRUE),
    sample(ct3_idx, n_doublets_13, replace = TRUE)
  )
)

doublet_counts <- lapply(1:nrow(doublet_pairs), function(i) {
  syn_data$counts[doublet_pairs[i, 1], ] +
    syn_data$counts[doublet_pairs[i, 2], ]
})
doublet_matrix <- do.call(rbind, doublet_counts)
all_counts <- rbind(syn_data$counts, doublet_matrix)
n_total <- nrow(all_counts)

# update obs
doublet_obs <- data.table::data.table(
  cell_id = sprintf("doublet_%04d", 1:n_doublets),
  cell_grp = "doublet",
  batch_index = 1,
  doublet = TRUE
)
new_obs <- data.table::rbindlist(list(
  syn_data$obs[, doublet := FALSE],
  doublet_obs
))

# test scrublet ----------------------------------------------------------------

## generate the object ---------------------------------------------------------

sc_object <- single_cell_exp(
  dir_data = tempdir()
)

# keep all cells for the sake of this
sc_object <- load_r_data(
  object = sc_object,
  counts = all_counts,
  obs = new_obs,
  var = syn_data$var,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 0L,
    min_lib_size = 0L,
    min_cells = 0L
  ),
  streaming = FALSE,
  .verbose = FALSE
)

## rust logic ------------------------------------------------------------------

test_wide_scrublet_params <- params_scrublet(
  log_transform = TRUE,
  mean_center = TRUE,
  normalise_variance = TRUE,
  no_pcs = 5L,
  expected_doublet_rate = 0.2,
  sim_doublet_ratio = 1.5,
  min_gene_var_pctl = 0.7,
  target_size = 1e4
)

scrublet_res = rs_sc_scrublet(
  f_path_gene = bixverse:::get_rust_count_gene_f_path(sc_object),
  f_path_cell = bixverse:::get_rust_count_cell_f_path(sc_object),
  cells_to_keep = get_cells_to_keep(sc_object),
  scrublet_params = test_wide_scrublet_params,
  seed = 42L,
  verbose = TRUE,
  streaming = FALSE
)

hist(scrublet_res$doublet_scores_obs)

hist(scrublet_res$doublet_scores_sim)

table(list(
  actual = new_obs$doublet,
  predicted = scrublet_res$predicted_doublets
))

metrics_helper(
  cm = table(new_obs$doublet, scrublet_res$predicted_doublets)
)

expect_true(
  current = checkmate::qtest(
    scrublet_res$predicted_doublets,
    sprintf("B%s", nrow(new_obs))
  ),
  info = "rust scrublet: predicted doublets correct return type"
)

expect_true(
  current = checkmate::qtest(
    scrublet_res$doublet_scores_obs,
    sprintf("N%s", nrow(new_obs))
  ),
  info = "rust scrublet: doublet_scores_obs correct return type"
)

# should be 1.5 times the length of the others...
expect_true(
  current = checkmate::qtest(
    scrublet_res$doublet_scores_sim,
    sprintf("N%s", nrow(new_obs) * test_wide_scrublet_params$sim_doublet_ratio)
  ),
  info = "rust scrublet: doublet_scores_sim correct return type"
)

expect_true(
  current = checkmate::qtest(
    scrublet_res$doublet_errors_obs,
    sprintf("N%s", nrow(new_obs))
  ),
  info = "rust scrublet: doublet_scores_sim correct return type"
)

expect_true(
  current = length(scrublet_res$predicted_doublets) == nrow(all_counts),
  info = "rust scrublet: correct dimensions"
)

metrics <- metrics_helper(
  cm = table(new_obs$doublet, scrublet_res$predicted_doublets)
)

expect_true(
  current = metrics["recall"] >= 0.7,
  info = "rust scrublet: 'good' recall on synthetic data"
)

expect_true(
  current = metrics["f1"] >= 0.7,
  info = "rust scrublet: 'good' recall on synthetic data"
)


## object ----------------------------------------------------------------------

obj_res <- scrublet_sc(
  sc_object,
  scrublet_params = test_wide_scrublet_params,
  .verbose = FALSE
)

hist(obj_res$doublet_scores_obs)

hist(obj_res$doublet_scores_sim)

hist(c(obj_res$doublet_scores_obs, obj_res$doublet_scores_sim))

plot(obj_res)

expect_true(
  current = checkmate::testClass(obj_res, "scrublet_res"),
  info = "S7 scrublet: the correct class is being returned"
)

expect_true(
  current = checkmate::qtest(attr(obj_res, "cell_indices"), "I+"),
  info = "S7 scrublet: the cell indices are stored on the class"
)

expect_equivalent(
  current = obj_res$predicted_doublets,
  target = scrublet_res$predicted_doublets,
  info = "S7 scrublet: no weird changes during generation (called doublets)"
)

expect_equivalent(
  current = obj_res$doublet_scores_obs,
  target = scrublet_res$doublet_scores_obs,
  info = "S7 scrublet: no weird changes during generation (obs scores)"
)

class(obj_res)
plot(obj_res)


# Scrublet simulation in R for debugging
scrublet_debug <- function(
  counts_matrix,
  n_doublets = 1000,
  target_size = 1e3,
  n_pcs = 30
) {
  # Step 1: Select HVG (using top 15% most variable)
  lib_sizes <- rowSums(counts_matrix)
  cpm <- sweep(counts_matrix, 1, lib_sizes, "/") * target_size
  gene_vars <- apply(cpm, 2, var)
  hvg_idx <- which(gene_vars > quantile(gene_vars, 0.85))

  cat("Selected", length(hvg_idx), "HVGs\n")

  # Step 2: Simulate doublets (add raw counts for HVG only)
  n_cells <- nrow(counts_matrix)
  doublet_pairs <- replicate(
    n_doublets,
    sample(n_cells, 2, replace = TRUE),
    simplify = FALSE
  )

  doublet_counts <- lapply(doublet_pairs, function(pair) {
    colSums(counts_matrix[pair, hvg_idx, drop = FALSE])
  })
  doublet_matrix <- do.call(rbind, doublet_counts)

  # Step 3: CPM normalize observed cells (HVG only)
  obs_hvg_counts <- counts_matrix[, hvg_idx]
  obs_lib_sizes <- rowSums(counts_matrix) # ALL genes
  obs_cpm <- sweep(obs_hvg_counts, 1, obs_lib_sizes, "/") * target_size

  # Step 4: CPM normalize doublets
  doublet_lib_sizes <- sapply(doublet_pairs, function(pair) {
    sum(lib_sizes[pair])
  })
  sim_cpm <- sweep(doublet_matrix, 1, doublet_lib_sizes, "/") * target_size

  # Step 5: Z-score using observed cell stats
  combined_cpm <- rbind(obs_cpm, sim_cpm)

  combined_cpm <- log1p(combined_cpm)

  gene_means <- colMeans(obs_cpm)
  gene_sds <- apply(obs_cpm, 2, sd)

  scaled <- sweep(sweep(combined_cpm, 2, gene_means, "-"), 2, gene_sds, "/")

  # Step 6: PCA
  pca_res <- prcomp(scaled, center = FALSE, scale. = FALSE, rank. = n_pcs)

  # Return results
  list(
    pca_scores = pca_res$x,
    n_obs = n_cells,
    n_sim = n_doublets,
    hvg_idx = hvg_idx,
    obs_cpm = obs_cpm,
    sim_cpm = sim_cpm,
    scaled = scaled
  )
}

# Run on your data
res <- scrublet_debug(
  as.matrix(all_counts),
  n_doublets = 1000,
  target_size = 1e4
)

# Visualize
library(ggplot2)
plot_data <- data.frame(
  PC1 = res$pca_scores[, 1],
  PC2 = res$pca_scores[, 2],
  type = c(rep("Observed", res$n_obs), rep("Simulated", res$n_sim)),
  is_doublet = c(new_obs$doublet[1:res$n_obs], rep(TRUE, res$n_sim))
)

ggplot(plot_data, aes(PC1, PC2, colour = type, shape = is_doublet)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(title = "Scrublet PCA: Observed vs Simulated Doublets")

# Check for NaN/Inf
cat("NaN in scaled:", sum(is.nan(res$scaled)), "\n")
cat("Inf in scaled:", sum(is.infinite(res$scaled)), "\n")
cat("Range of PCA scores:", range(res$pca_scores), "\n")

head(plot_data, 15)

tail(plot_data, 15)


# Check intermediate steps
cat("=== Data Quality Checks ===\n")
cat("Observed cells sparsity:", mean(res$obs_cpm == 0), "\n")
cat("Simulated cells sparsity:", mean(res$sim_cpm == 0), "\n")
cat("Observed CPM range:", range(res$obs_cpm), "\n")
cat("Simulated CPM range:", range(res$sim_cpm), "\n")

# Check if scaling introduced issues
cat("\n=== Scaling Stats ===\n")
cat("Gene means range:", range(colMeans(res$obs_cpm)), "\n")
cat("Gene SDs range:", range(apply(res$obs_cpm, 2, sd)), "\n")
cat(
  "Genes with SD < 0.01:",
  sum(apply(res$obs_cpm, 2, sd) < 0.01),
  "/",
  ncol(res$obs_cpm),
  "\n"
)

# Check the original count matrix
cat("\n=== Original Data ===\n")
cat("Count matrix sparsity:", mean(as.matrix(all_counts) == 0), "\n")
cat("Lib sizes range:", range(rowSums(as.matrix(all_counts))), "\n")
cat("Mean counts per gene:", mean(colSums(as.matrix(all_counts))), "\n")

# Visualize first few PCs
pairs(
  res$pca_scores[, 1:3],
  col = ifelse(plot_data$type == "Observed", "blue", "red"),
  pch = 16,
  cex = 0.5,
  main = "First 3 PCs"
)

# Check variance explained
pca_full <- prcomp(res$scaled, center = FALSE, scale. = FALSE)
plot(
  summary(pca_full)$importance[2, 1:10],
  type = "b",
  main = "Variance Explained by PC",
  xlab = "PC",
  ylab = "Proportion of Variance"
)

# Most importantly: Check if your synthetic data has cell type structure
if (exists("new_obs")) {
  obs_only <- plot_data[plot_data$type == "Observed", ]
  obs_only$cell_type <- new_obs$cell_grp[1:nrow(obs_only)]

  ggplot(obs_only, aes(PC1, PC2, colour = cell_type)) +
    geom_point(alpha = 0.6) +
    theme_bw() +
    labs(title = "Do cell types separate in PCA?")
}


# check some other data --------------------------------------------------------

demuxlet_result <- fread("~/Downloads/demuxlet_PBMCs/demuxlet_calls.tsv")

sc_object_pmbc <- single_cell_exp(
  dir_data = tempdir()
)

mtx_params <- params_sc_mtx_io(
  path_mtx = path.expand("~/Downloads/demuxlet_PBMCs/output.mtx"),
  path_obs = path.expand("~/Downloads/demuxlet_PBMCs/GSM2560248_barcodes.tsv"),
  path_var = path.expand("~/Downloads/demuxlet_PBMCs/GSM2560248_genes.tsv"),
  cells_as_rows = TRUE,
  has_hdr = FALSE
)

sc_object_pmbc <- load_mtx(
  sc_object_pmbc,
  sc_mtx_io_param = mtx_params,
  sc_qc_param = params_sc_min_quality()
)

sc_object_pmbc@dims

devtools::load_all()

rextendr::document()

scrublet_params = params_scrublet(
  log_transform = FALSE,
  mean_center = TRUE,
  normalise_variance = TRUE,
  expected_doublet_rate = 0.05,
  dist_metric = "euclidean",
  target_size = 1e4,
  sim_doublet_ratio = 2.0
)

scrublet_res = rs_sc_scrublet(
  f_path_gene = bixverse:::get_rust_count_gene_f_path(sc_object_pmbc),
  f_path_cell = bixverse:::get_rust_count_cell_f_path(sc_object_pmbc),
  cells_to_keep = get_cells_to_keep(sc_object_pmbc),
  scrublet_params = scrublet_params,
  seed = 1210103L,
  verbose = TRUE,
  streaming = FALSE
)

scrublet_res$doublet_scores_obs[1:10]

pmbc_scrublet <- scrublet_sc(
  sc_object_pmbc,
  scrublet_params = scrublet_params,
  .verbose = TRUE
)

plot(pmbc_scrublet)

pmbc_scrublet <- call_doublets_manual(pmbc_scrublet, threshold = 0.05)

pmbc_scrublet$doublet_scores_obs[1:10]

plot(pmbc_scrublet)

obs <- sc_object_pmbc[[]]
obs[, doublet := pmbc_scrublet$predicted_doublets]

demuxlet_result_combined <- merge(
  obs,
  demuxlet_result,
  by.x = "cell_id",
  by.y = "Barcode"
)[Call != "AMB"][, doublet_demuxlet := Call == "DBL"]

table(
  list(
    predicted = demuxlet_result_combined$doublet,
    actual = demuxlet_result_combined$doublet_demuxlet
  )
)

metrics_helper(table(
  list(
    predicted = demuxlet_result_combined$doublet,
    actual = demuxlet_result_combined$doublet_demuxlet
  )
))


sum(demuxlet_result_combined$doublet_demuxlet) /
  nrow(demuxlet_result_combined)

print(table(
  demuxlet_result_combined$doublet,
  demuxlet_result_combined$Call
))
