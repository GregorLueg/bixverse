# sc processing ----------------------------------------------------------------

set.seed(123L)

test_temp_dir <- file.path(
  tempdir(),
  "processing"
)

dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## testing parameters ----------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
# hvg
hvg_to_keep <- 30L
# pca
no_pcs <- 15L

## synthetic test data ---------------------------------------------------------

single_cell_test_data <- generate_single_cell_test_data()

genes_pass <- which(
  Matrix::colSums(single_cell_test_data$counts != 0) >= min_cells_exp
)

cells_pass <- which(
  (Matrix::rowSums(single_cell_test_data$counts[, genes_pass]) >=
    min_lib_size) &
    (Matrix::rowSums(single_cell_test_data$counts[, genes_pass] != 0) >=
      min_genes_exp)
)

expect_true(
  current = length(genes_pass) > 80 & length(genes_pass) != 100,
  info = "sc processing - sensible amount of genes pass"
)

expect_true(
  current = length(cells_pass) > 800 & length(cells_pass) != 1000,
  info = "sc processing - sensible amount of cells pass"
)

counts_filtered <- single_cell_test_data$counts[cells_pass, genes_pass]

sc_qc_param <- params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

## underlying class ------------------------------------------------------------

sc_object <- SingleCells(dir_data = test_temp_dir)

sc_object <- load_r_data(
  object = sc_object,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  sc_qc_param = sc_qc_param,
  streaming = 0L,
  .verbose = FALSE
)

# tests ------------------------------------------------------------------------

## test filter column and cells to get -----------------------------------------

expect_true(
  current = checkmate::qtest(unlist(sc_object[["to_keep"]]), "B+"),
  info = "to_keep column added"
)

expect_true(
  current = checkmate::qtest(get_cells_to_keep(sc_object), "I+"),
  info = "cells_to_keep loaded independently of setting"
)

## function warnings -----------------------------------------------------------

expect_warning(
  current = get_hvg(sc_object),
  info = "warning that hvgs can be found"
)

expect_warning(
  current = calculate_pca_sc(sc_object, no_pcs = 10L),
  info = "warning that no HVGs are detected"
)

expect_warning(
  current = find_neighbours_sc(object = sc_object),
  info = "warning that no PCA data are detected"
)

expect_warning(
  current = find_clusters_sc(sc_object),
  info = "warning that no kNN/sNN data was found"
)

## gene set proportions --------------------------------------------------------

### top proportion n genes -----------------------------------------------------

top_n <- c(5L, 10L, 25L)

proportions_in_r <- sapply(top_n, function(x) {
  apply(as.matrix(counts_filtered), 1, function(row) {
    sum(sort(as.matrix(row), decreasing = TRUE)[1:x])
  }) /
    rowSums(as.matrix(counts_filtered))
})

sc_object <- top_genes_perc_sc(sc_object, top_n_vals = top_n, .verbose = FALSE)

expect_equivalent(
  current = unlist(sc_object[["top_5_genes_percentage"]]),
  target = proportions_in_r[, 1],
  tolerance = 1e-7,
)

expect_equivalent(
  current = unlist(sc_object[["top_10_genes_percentage"]]),
  target = proportions_in_r[, 2],
  tolerance = 1e-7,
)

expect_equivalent(
  current = unlist(sc_object[["top_25_genes_percentage"]]),
  target = proportions_in_r[, 3],
  tolerance = 1e-7,
)

#### streaming -----------------------------------------------------------------

sc_object <- top_genes_perc_sc(
  sc_object,
  top_n_vals = top_n,
  streaming = TRUE,
  .verbose = FALSE
)

expect_equivalent(
  current = unlist(sc_object[["top_5_genes_percentage"]]),
  target = proportions_in_r[, 1],
  tolerance = 1e-7,
)

expect_equivalent(
  current = unlist(sc_object[["top_10_genes_percentage"]]),
  target = proportions_in_r[, 2],
  tolerance = 1e-7,
)

expect_equivalent(
  current = unlist(sc_object[["top_25_genes_percentage"]]),
  target = proportions_in_r[, 3],
  tolerance = 1e-7,
)

### gene set -------------------------------------------------------------------

gs_of_interest <- list(
  gs_1 = c("gene_001", "gene_002", "gene_003", "gene_004"),
  gs_2 = c("gene_096", "gene_097", "gene_100")
)

total_size <- Matrix::rowSums(counts_filtered)

props_gs_1 <- Matrix::rowSums(counts_filtered[, gs_of_interest$gs_1]) /
  total_size
props_gs_2 <- Matrix::rowSums(counts_filtered[, gs_of_interest$gs_2]) /
  total_size

sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  .verbose = FALSE
)

expect_equivalent(
  current = unlist(sc_object[["gs_1"]]),
  target = props_gs_1,
  tolerance = 10e-7,
  info = "gene proportion calculations gene set 1"
)

expect_equivalent(
  current = unlist(sc_object[["gs_2"]]),
  target = props_gs_2,
  tolerance = 10e-7,
  info = "gene proportion calculations gene set 2"
)

sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  .verbose = FALSE
)

expect_true(
  current = all(!c("gs_1.1", "gs_2.1") %in% colnames(sc_object[[]])),
  info = "overwriting of obs data works"
)

#### streaming -----------------------------------------------------------------

sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  streaming = TRUE,
  .verbose = FALSE
)

expect_equivalent(
  current = unlist(sc_object[["gs_1"]]),
  target = props_gs_1,
  tolerance = 10e-7,
  info = "gene proportion calculations gene set 1 - streaming version"
)

expect_equivalent(
  current = unlist(sc_object[["gs_2"]]),
  target = props_gs_2,
  tolerance = 10e-7,
  info = "gene proportion calculations gene set 2 - streaming version"
)

## cells to keep logic ---------------------------------------------------------

threshold <- 0.05

cells_to_keep <- sc_object[[]][gs_2 < threshold, cell_id]

expect_true(
  current = length(cells_to_keep) > 600,
  info = "sensible cell filtering based on the threshold"
)

sc_object <- set_cells_to_keep(sc_object, cells_to_keep)

expect_true(
  current = all(
    unlist(sc_object[["cell_id"]]) == cells_to_keep
  ),
  info = "setting cells to keep removes them from the obs table"
)

cell_names_filtered <- get_cell_names(sc_object, filtered = TRUE)

expect_true(
  current = all(cell_names_filtered == cells_to_keep),
  info = "filter flag on cell names works"
)

counts_more_filtered <- counts_filtered[
  which(props_gs_2 < threshold),
]

expect_equivalent(
  current = get_cells_to_keep(sc_object),
  target = which(props_gs_2 < threshold) - 1, # zero indexed in Rust
  info = "cell to keep logic"
)

# the logic of retrieving cells will become more complicated here...

expect_equal(
  current = sc_object[,, use_cells_to_keep = TRUE],
  target = counts_more_filtered,
  info = "counts after setting cells to keep and using the parameter"
)

expect_equal(
  current = sc_object[,, use_cells_to_keep = FALSE],
  target = counts_filtered,
  info = "counts after setting cells to keep and NOT using the parameter"
)

## check addition of new data --------------------------------------------------

new_data <- rep("random_new_data", length(cells_to_keep))

new_data_list <- list(
  "other_random_data" = rep("A", length(cells_to_keep)),
  "even_different_random_data" = seq_len(length(cells_to_keep))
)

sc_object[["random_new_data"]] <- new_data

sc_object[[names(new_data_list)]] <- new_data_list

sc_object[[c("new_name_a", "new_name_b")]] <- new_data_list

expect_true(
  current = unique(unlist(sc_object[[c("random_new_data")]])) ==
    "random_new_data",
  info = "obs table addition worked - single value"
)

expect_true(
  current = unique(unlist(sc_object[[c("other_random_data")]])) == "A",
  info = "obs table addition worked - from list string"
)

expect_true(
  current = all(
    unlist(sc_object[[c("even_different_random_data")]]) ==
      seq_len(length(cells_to_keep))
  ),
  info = "obs table addition worked - from list numeric"
)

expect_true(
  current = unique(unlist(sc_object[[c("new_name_a")]])) == "A",
  info = "obs table addition worked - from list string - renamed"
)

expect_true(
  current = all(
    unlist(sc_object[[c("new_name_b")]]) == seq_len(length(cells_to_keep))
  ),
  info = "obs table addition worked - from list string - renamed"
)

# error handling

new_data_bad <- c("A")

new_data_bad_v2 <- list(
  "x" = c(1, 2, 3),
  "y" = letters
)

expect_error(
  current = {
    sc_object[["random_new_data"]] <- new_data_bad
  },
  info = "addition of wrong length throws the right error"
)

expect_error(
  current = {
    sc_object[[names(new_data_bad_v2)]] <- new_data_bad_v2
  },
  info = "addition of list with random lengths throws an error"
)

### test that the original data is still in there ------------------------------

expect_true(
  current = nrow(get_sc_obs(sc_object)) == sc_object@dims[1],
  info = "original rows are still in the DB"
)

## hvg selection ---------------------------------------------------------------

### vst version ----------------------------------------------------------------

#### r version -----------------------------------------------------------------

n_cells <- nrow(counts_more_filtered)

# calculate mean and variance
mean_values_r <- Matrix::colMeans(counts_more_filtered)
var_values_r <- matrixStats::colVars(as.matrix(counts_more_filtered)) *
  ((n_cells - 1) /
    n_cells)

valid_genes <- var_values_r > 0 & mean_values_r > 0

mean_log10 <- log10(mean_values_r[valid_genes])
var_log10 <- log10(var_values_r[valid_genes])

loess_fit <- loess(var_log10 ~ mean_log10, span = 0.3, degree = 2)
var_expected_log10 <- predict(loess_fit, mean_log10)
var_expected <- 10^var_expected_log10

clip_max <- sqrt(n_cells)

var_std <- sapply(1:ncol(counts_more_filtered), function(i) {
  gene_vals <- counts_more_filtered[, i]
  mean_i <- mean_values_r[i]
  sd_expected <- sqrt(var_expected[i])

  standardised <- (gene_vals - mean_i) / sd_expected

  standardised <- pmin(pmax(standardised, -clip_max), clip_max)

  # Variance of standardized values
  # population variance
  var(standardised) * ((n_cells - 1) / n_cells)
})

var_std[is.na(var_std)] <- 0

# need to 0-index
hvg_r <- order(var_std, decreasing = TRUE)[1:hvg_to_keep] - 1

#### rust part -----------------------------------------------------------------

##### direct load --------------------------------------------------------------

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

var_data <- get_sc_var(sc_object)

expect_equivalent(
  current = var_data$mean,
  target = mean_values_r,
  tolerance = 10e-7,
  info = "Correct mean calculations for genes"
)

expect_equivalent(
  current = var_data$var,
  target = var_values_r,
  tolerance = 10e-7,
  info = "Correct variance calculations for genes"
)

expect_equivalent(
  current = var_data$var_std,
  target = var_std,
  tolerance = 10e-3,
  info = "Correct standardised variance calculations for genes"
)

hvg_rs <- get_hvg(sc_object)

expect_true(
  current = length(intersect(hvg_r, hvg_rs)) == hvg_to_keep,
  info = "Overlap in the detected HVGs"
)

##### streaming version --------------------------------------------------------

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  streaming = TRUE,
  .verbose = FALSE
)

var_data <- get_sc_var(sc_object)

expect_equivalent(
  current = var_data$mean,
  target = mean_values_r,
  tolerance = 10e-7,
  info = "Correct mean calculations for genes"
)

expect_equivalent(
  current = var_data$var,
  target = var_values_r,
  tolerance = 10e-7,
  info = "Correct variance calculations for genes"
)

expect_equivalent(
  current = var_data$var_std,
  target = var_std,
  tolerance = 10e-3,
  info = "Correct standardised variance calculations for genes"
)

hvg_rs <- get_hvg(sc_object)

expect_true(
  current = length(intersect(hvg_r, hvg_rs)) == hvg_to_keep,
  info = "Overlap in the detected HVGs"
)

### dispersion version --------------------------------------------------------

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  hvg_params = params_sc_hvg(method = "dispersion"),
  .verbose = FALSE
)

var_data_disp <- get_sc_var(sc_object)

expect_true(
  current = all(
    c("mean", "dispersion", "dispersion_scaled", "bin") %in%
      names(var_data_disp)
  ),
  info = "sc hvg dispersion - var table populated with dispersion metrics"
)

hvg_disp <- get_hvg(sc_object)

expect_true(
  current = checkmate::qtest(hvg_disp, sprintf("I%i", hvg_to_keep)),
  info = "sc hvg dispersion - correct number of HVGs returned"
)

# signal genes (1-30) should dominate non-signal genes (31-100)
expect_true(
  current = length(intersect(hvg_disp + 1, 1:30)) >=
    length(intersect(hvg_disp + 1, 31:100)),
  info = "sc hvg dispersion - signal genes dominate the HVGs"
)

#### streaming version --------------------------------------------------------

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  hvg_params = params_sc_hvg(method = "dispersion"),
  streaming = TRUE,
  .verbose = FALSE
)

var_data_disp_stream <- get_sc_var(sc_object)

expect_equivalent(
  current = var_data_disp_stream$mean,
  target = var_data_disp$mean,
  tolerance = 1e-6,
  info = "sc hvg dispersion - streaming matches non-streaming (mean)"
)

expect_equivalent(
  current = var_data_disp_stream$dispersion,
  target = var_data_disp$dispersion,
  tolerance = 1e-5,
  info = "sc hvg dispersion - streaming matches non-streaming (dispersion)"
)

expect_equivalent(
  current = var_data_disp_stream$dispersion_scaled,
  target = var_data_disp$dispersion_scaled,
  tolerance = 1e-5,
  info = "sc hvg dispersion - streaming matches non-streaming (scaled)"
)

expect_true(
  current = length(intersect(get_hvg(sc_object), hvg_disp)) == hvg_to_keep,
  info = "sc hvg dispersion - streaming returns the same HVGs"
)

### meanvarbin version -------------------------------------------------------

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  hvg_params = params_sc_hvg(
    method = "meanvarbin",
    num_bin = 5L
  ),
  .verbose = FALSE
)

var_data_mvb <- get_sc_var(sc_object)

expect_true(
  current = all(
    c("mean", "dispersion", "dispersion_scaled", "bin") %in% names(var_data_mvb)
  ),
  info = "sc hvg meanvarbin - var table populated with metrics"
)

hvg_mvb <- get_hvg(sc_object)

expect_true(
  current = checkmate::qtest(hvg_mvb, sprintf("I%i", hvg_to_keep)),
  info = "sc hvg meanvarbin - correct number of HVGs returned"
)

# this one behaves worse than the others due to the weirdness in the synthetic
# data
expect_true(
  current = length(intersect(hvg_mvb + 1, 1:30)) >= 10L,
  info = "sc hvg meanvarbin - signal genes dominate the HVGs"
)

#### streaming version --------------------------------------------------------

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  hvg_params = params_sc_hvg(method = "meanvarbin", num_bin = 5L),
  streaming = TRUE,
  .verbose = FALSE
)

var_data_mvb_stream <- get_sc_var(sc_object)

expect_equivalent(
  current = var_data_mvb_stream$mean,
  target = var_data_mvb$mean,
  tolerance = 1e-6,
  info = "sc hvg meanvarbin - streaming matches non-streaming (mean)"
)

expect_equivalent(
  current = var_data_mvb_stream$dispersion_scaled,
  target = var_data_mvb$dispersion_scaled,
  tolerance = 1e-5,
  info = "sc hvg meanvarbin - streaming matches non-streaming (scaled)"
)

expect_true(
  current = length(intersect(get_hvg(sc_object), hvg_mvb)) == hvg_to_keep,
  info = "sc hvg meanvarbin - streaming returns the same HVGs"
)

### get_hvg_data_sc -----------------------------------------------------------

# reset to vst so downstream PCA tests remain unchanged
sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

current_hvg <- get_hvg(sc_object)
current_var <- get_sc_var(sc_object)

#### full cells_to_keep -------------------------------------------------------

hvg_dt_full <- get_hvg_data_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testDataTable(hvg_dt_full),
  info = "get_hvg_data_sc - returns a data.table"
)

expect_true(
  current = all(
    c(
      "gene_idx",
      "gene_id",
      "mean",
      "var",
      "var_std",
      "is_hvg",
      "hvg_rank"
    ) %in%
      names(hvg_dt_full)
  ),
  info = "get_hvg_data_sc - expected columns present"
)

expect_true(
  current = sum(hvg_dt_full$is_hvg) == hvg_to_keep,
  info = "get_hvg_data_sc - correct number of HVGs flagged"
)

expect_true(
  current = all(!is.na(hvg_dt_full[is_hvg == TRUE, hvg_rank])) &
    all(is.na(hvg_dt_full[is_hvg == FALSE, hvg_rank])),
  info = "get_hvg_data_sc - hvg_rank NA exactly when is_hvg is FALSE"
)

# cell_ids = NULL should give HVGs equivalent to those stored on the object
expect_true(
  current = length(intersect(
    hvg_dt_full[is_hvg == TRUE, gene_idx],
    current_hvg + 1
  )) ==
    hvg_to_keep,
  info = "get_hvg_data_sc - NULL cell_ids matches state-mutating version"
)

#### biological subset -------------------------------------------------------

ct1_cells <- sc_object[[]][cell_grp == "cell_type_1", cell_id]

expect_true(
  current = length(ct1_cells) > 50,
  info = "get_hvg_data_sc - sensible number of cell_type_1 cells for subset test"
)

hvg_dt_subset <- get_hvg_data_sc(
  object = sc_object,
  cell_ids = ct1_cells,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

expect_true(
  current = sum(hvg_dt_subset$is_hvg) == hvg_to_keep,
  info = "get_hvg_data_sc - subset returns correct number of HVGs"
)

# means/vars computed on the subset must differ from the full data
expect_false(
  current = isTRUE(all.equal(hvg_dt_subset$mean, hvg_dt_full$mean)),
  info = "get_hvg_data_sc - subset stats differ from full stats"
)

#### state preservation -----------------------------------------------------

expect_equal(
  current = get_hvg(sc_object),
  target = current_hvg,
  info = "get_hvg_data_sc - object HVGs unchanged after call"
)

expect_equivalent(
  current = get_sc_var(sc_object)$var_std,
  target = current_var$var_std,
  info = "get_hvg_data_sc - object var table unchanged after call"
)

## pca -------------------------------------------------------------------------

### seurat version -------------------------------------------------------------

#### r -------------------------------------------------------------------------

pca_input <- as.matrix(sc_object[,
  as.integer(get_hvg(sc_object) + 1),
  assay = "norm",
  return_format = "gene",
  use_cells_to_keep = TRUE
])

scaled_data <- scale(pca_input)

pca_r <- prcomp(pca_input, scale. = TRUE)

expected_names <- get_gene_names(sc_object)[hvg_r + 1]
actual_names <- colnames(pca_input)

#### rust ----------------------------------------------------------------------

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = no_pcs,
  pca_params = params_sc_pca(
    mean_center = TRUE,
    normalise_variance = TRUE,
    clr = FALSE,
    randomised = FALSE
  ),
  .verbose = FALSE
)

expect_true(
  current = all.equal(
    abs(diag(cor(get_pca_factors(sc_object)[, 1:no_pcs], pca_r$x[, 1:no_pcs]))),
    rep(1, no_pcs),
    tolerance = 1e-8
  ),
  info = "PCA on single cell data compared to R"
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = no_pcs,
  pca_params = params_sc_pca(
    mean_center = TRUE,
    normalise_variance = TRUE,
    clr = FALSE,
    randomised = TRUE
  ),
  .verbose = FALSE
)

expect_true(
  current = all.equal(
    abs(diag(cor(get_pca_factors(sc_object)[, 1:no_pcs], pca_r$x[, 1:no_pcs]))),
    rep(1, no_pcs),
    tolerance = 1e-8
  ),
  info = "PCA on single cell data compared to R (randomised SVD)"
)

#### scaling within the rust function ------------------------------------------

# test if the scaling results in the same data
zeallot::`%<-%`(
  c(pca_factors, pca_loadings, pca_eigenvals, scaled),
  rs_sc_pca(
    f_path_gene = bixverse:::get_rust_count_gene_f_path(sc_object),
    f_path_cell = bixverse:::get_rust_count_gene_f_path(sc_object),
    no_pcs = no_pcs,
    pca_params = params_sc_pca(
      mean_center = TRUE,
      normalise_variance = TRUE,
      clr = FALSE,
      randomised = FALSE
    ),
    cell_indices = get_cells_to_keep(sc_object),
    gene_indices = get_hvg(sc_object),
    seed = 42L,
    return_scaled = TRUE,
    verbose = 0L
  )
)

expect_equivalent(
  current = scaled,
  target = scaled_data,
  tolerance = 1e-6,
  info = "scaling behaves"
)

### PFlogPF version ------------------------------------------------------------

# CLR-type normalisation

#### r -------------------------------------------------------------------------

# pre-calculated everything in R

raw_counts <- as.matrix(sc_object[,,
  assay = "raw",
  return_format = "cell",
  use_cells_to_keep = TRUE
])

lib_sizes <- rowSums(raw_counts)
u <- sweep(raw_counts, 1, lib_sizes, "/")
log1p_u <- log1p(u)

clr_offsets <- rowMeans(log1p_u)

clr_matrix <- sweep(log1p_u, 1, clr_offsets, "-")
pca_input_clr <- clr_matrix[, as.integer(get_hvg(sc_object) + 1)]

pca_r_clr <- prcomp(pca_input_clr, center = FALSE, scale. = FALSE)

pca_r_clr_scaled <- prcomp(pca_input_clr, center = TRUE, scale. = TRUE)

#### rust ----------------------------------------------------------------------

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = no_pcs,
  pca_params = params_sc_pca(
    mean_center = FALSE,
    normalise_variance = FALSE,
    clr = TRUE,
    randomised = FALSE,
    size_factor = 1e3
  ),
  .verbose = FALSE
)

expect_true(
  current = all.equal(
    abs(diag(cor(
      get_pca_factors(sc_object)[, 1:no_pcs],
      pca_r_clr$x[, 1:no_pcs]
    ))),
    rep(1, no_pcs),
    tolerance = 1e-5
  ),
  info = "PFlogPF PCA on single cell data compared to R"
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = no_pcs,
  pca_params = params_sc_pca(
    mean_center = TRUE,
    normalise_variance = TRUE,
    clr = TRUE,
    randomised = FALSE,
    size_factor = 1e3
  ),
  .verbose = FALSE
)

expect_true(
  current = all.equal(
    abs(diag(cor(
      get_pca_factors(sc_object)[, 1:no_pcs],
      pca_r_clr_scaled$x[, 1:no_pcs]
    ))),
    rep(1, no_pcs),
    tolerance = 1e-5
  ),
  info = "PFlogPF PCA (scaled) on single cell data compared to R"
)

#### compare scaled data -------------------------------------------------------

# test if the scaling results in the same data
zeallot::`%<-%`(
  c(pca_factors_clr, pca_loadings_clr, pca_eigenvals_clr, scaled_clr),
  rs_sc_pca(
    f_path_gene = bixverse:::get_rust_count_gene_f_path(sc_object),
    f_path_cell = bixverse:::get_rust_count_cell_f_path(sc_object),
    no_pcs = no_pcs,
    pca_params = params_sc_pca(
      mean_center = FALSE,
      normalise_variance = FALSE,
      clr = TRUE,
      randomised = FALSE,
      size_factor = 1e3
    ),
    cell_indices = get_cells_to_keep(sc_object),
    gene_indices = get_hvg(sc_object),
    seed = 42L,
    return_scaled = TRUE,
    verbose = 0L
  )
)

expect_equivalent(
  current = scaled_clr,
  target = pca_input_clr,
  tolerance = 1e-3, # FP16
  info = "clr transformation behaves"
)

### rerun pca with normal parameters -------------------------------------------

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = no_pcs,
  pca_params = params_sc_pca(),
  .verbose = FALSE
)

## knn and snn -----------------------------------------------------------------

sc_object <- find_neighbours_sc(
  sc_object,
  .verbose = FALSE
)

expect_equal(
  current = dim(get_knn_mat(sc_object)),
  target = c(682, 15),
  info = "kNN matrix correctly returned"
)

expect_true(
  current = class(get_snn_graph(sc_object)) == "igraph",
  info = "igraph correctly returned"
)

## community detection ---------------------------------------------------------

sc_object <- find_clusters_sc(sc_object)

cell_grps <- unlist(sc_object[["cell_grp"]])
leiden_clusters <- unlist(sc_object[["leiden_clustering"]])

cell_grps <- unlist(sc_object[["cell_grp"]])

f1_scores <- f1_score_confusion_mat(cell_grps, leiden_clusters)

expect_true(
  current = all(f1_scores > 0.95),
  info = "leiden clustering identifies the cell groups"
)

## fast clustering -------------------------------------------------------------

### simple version -------------------------------------------------------------

fast_cluster_res <- fast_cluster_sc(
  object = sc_object,
  .verbose = FALSE
)

obs_fc <- get_data(fast_cluster_res)

expect_true(
  current = checkmate::testDataTable(obs_fc),
  info = "fast clustering: data.table returned"
)

expect_warning(
  current = get_centroids(fast_cluster_res),
  info = "fast clustering: warning without k-means data"
)

expect_warning(
  current = get_kmeans_clusters(fast_cluster_res),
  info = "fast clustering: warning without k-means data (2)"
)

### with k-means clusters ------------------------------------------------------

fast_cluster_res <- fast_cluster_sc(
  object = sc_object,
  return_kmeans = TRUE,
  n_centroids = 30L,
  .verbose = FALSE
)

centroids <- get_centroids(fast_cluster_res)

kmeans_clusters <- get_kmeans_clusters(fast_cluster_res)

expect_true(
  current = checkmate::testMatrix(centroids, nrow = 30L, ncols = no_pcs),
  info = "fast clustering: centroids returned"
)

expect_true(
  current = checkmate::qtest(kmeans_clusters, "I+"),
  info = "fast clustering: k-means clusters returned"
)

### grid version ---------------------------------------------------------------

### check the DB structure -----------------------------------------------------

expect_true(
  current = all(is.na(get_sc_obs(sc_object)[!(to_keep), leiden_clustering])),
  info = "leiden clustering for cells that are not kept should be NA"
)

## dges ------------------------------------------------------------------------

### between two groups ---------------------------------------------------------

cell_names_1 <- sc_object[[]][cell_grp == "cell_type_1", cell_id]
cell_names_2 <- sc_object[[]][cell_grp == "cell_type_2", cell_id]

expect_error(
  current = find_markers_sc(
    object = sc_object,
    cells_1 = c(cell_names_1, "x"),
    cells_2 = cell_names_2
  ),
  info = "error if weird cells are selected"
)

dge_test <- find_markers_sc(
  object = sc_object,
  cells_1 = cell_names_1,
  cells_2 = cell_names_2,
  .verbose = FALSE
)

expected_upregulated <- sprintf("gene_%03d", 1:10)
expected_downregulated <- sprintf("gene_%03d", 11:20)

expect_true(
  current = all(
    expected_upregulated %in% dge_test[lfc > 0 & fdr <= 0.05, gene_id]
  ),
  info = "all expected up-regulated genes identified"
)

expect_true(
  current = all(
    expected_downregulated %in% dge_test[lfc < 0 & fdr <= 0.05, gene_id]
  ),
  info = "all expected down-regulated genes identified"
)

### find all markers -----------------------------------------------------------

dge_test_2 <- find_all_markers_sc(
  object = sc_object,
  column_of_interest = "leiden_clustering",
  .verbose = FALSE
)

all_cell_markers <- sprintf("gene_%03d", 1:30)

expect_true(
  current = all(
    all_cell_markers %in% dge_test_2[fdr <= 0.05, gene_id]
  ),
  info = "all expected cell markers identified"
)

## sparse pca ------------------------------------------------------------------

### direct rust functions ------------------------------------------------------

zeallot::`%<-%`(
  c(sparse_pca_factors, sparse_pca_loadings, sparse_pca_eigenvals),
  rs_sc_pca_sparse(
    f_path_gene = bixverse:::get_rust_count_gene_f_path(sc_object),
    f_path_cell = bixverse:::get_rust_count_cell_f_path(sc_object),
    no_pcs = no_pcs,
    pca_params = params_sc_pca(randomised = FALSE),
    cell_indices = get_cells_to_keep(sc_object),
    gene_indices = get_hvg(sc_object),
    seed = 42L,
    verbose = 0L
  )
)

zeallot::`%<-%`(
  c(sparse_pca_factors_rnd, sparse_pca_loadings_rnd, sparse_pca_eigenvals_rnd),
  rs_sc_pca_sparse(
    f_path_gene = bixverse:::get_rust_count_gene_f_path(sc_object),
    f_path_cell = bixverse:::get_rust_count_cell_f_path(sc_object),
    no_pcs = no_pcs,
    pca_params = params_sc_pca(randomised = FALSE),
    cell_indices = get_cells_to_keep(sc_object),
    gene_indices = get_hvg(sc_object),
    seed = 42L,
    verbose = 0L
  )
)

expect_equal(
  current = abs(diag(cor(sparse_pca_factors, sparse_pca_factors_rnd))),
  target = rep(1, 15),
  tolerance = 1e-7,
  into = "sparse svd - randomised and normal return very similar results"
)

expect_equal(
  current = abs(diag(cor(pca_factors, sparse_pca_factors))),
  target = rep(1, no_pcs),
  info = "sparse PCA returns the same as dense PCA"
)

### object method --------------------------------------------------------------

# more PCs needed for sparse SVD to recover all of the data set structure

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = no_pcs,
  pca_params = params_sc_pca(randomised = FALSE),
  sparse_svd = TRUE,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testMatrix(
    x = get_pca_factors(sc_object),
    mode = "numeric",
    row.names = "named",
    col.names = "named",
    ncol = no_pcs
  ),
  info = "sparse PCA correctly added to object and returned"
)

### test synthetic data set structure recovery ---------------------------------

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(),
  .verbose = FALSE
)

sc_object <- find_clusters_sc(sc_object)

cell_grps <- unlist(sc_object[["cell_grp"]])
leiden_clusters <- unlist(sc_object[["leiden_clustering"]])

cell_grps <- unlist(sc_object[["cell_grp"]])

f1_scores <- f1_score_confusion_mat(cell_grps, leiden_clusters)

expect_true(
  current = all(f1_scores > 0.95),
  info = "leiden clustering on top of sparse SVD identifies the cell groups"
)

# clean up ---------------------------------------------------------------------

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
