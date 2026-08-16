# sc aggregations --------------------------------------------------------------

source("helper_sc.R", local = TRUE)

library(magrittr)

test_temp_dir <- sc_test_dir("sc_aggregations")

## testing parameters ----------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
# hvg
hvg_to_keep <- 30L
# pca
no_pcs <- 10L

## synthetic test data ---------------------------------------------------------

single_cell_test_data <- sc_test_fixture(
  min_lib_size = min_lib_size,
  min_genes_exp = min_genes_exp,
  min_cells_exp = min_cells_exp,
  hvg_to_keep = hvg_to_keep,
  no_pcs = no_pcs
)

genes_pass <- single_cell_test_data$genes_pass
cells_pass <- single_cell_test_data$cells_pass

sc_qc_param <- sc_test_qc_params(single_cell_test_data, target_size = 1000)

## underlying class ------------------------------------------------------------

sc_object <- sc_test_prepped(
  sc_test_object(test_temp_dir, single_cell_test_data),
  single_cell_test_data
)

## pseudo bulking --------------------------------------------------------------

pseudobulk_list <- sc_object[[c("cell_grp", "cell_id")]] %$%
  split(cell_id, cell_grp)

pseudobulk_list_weird <- pseudobulk_list
pseudobulk_list_weird[[3]] <- c(pseudobulk_list_weird[[3]], "x")

### warnings / errors ----------------------------------------------------------

expect_error(
  current = suppressWarnings(get_pseudobulked_sc(
    object = sc_object,
    cell_list = list("a" = "a", "b" = "b"),
    return_format = "sparse",
    assay = "norm"
  )),
  info = "expect an error if there are zero matching cell indices"
)

expect_warning(
  current = get_pseudobulked_sc(
    object = sc_object,
    cell_list = pseudobulk_list_weird,
    return_format = "sparse",
    assay = "norm"
  ),
  info = "expect a warning if some identifiers cannot be matched"
)

### pseudo-bulks ---------------------------------------------------------------

matrix_raw <- get_pseudobulked_sc(
  object = sc_object,
  cell_list = pseudobulk_list
)

matrix_norm <- get_pseudobulked_sc(
  object = sc_object,
  cell_list = pseudobulk_list,
  assay = "norm"
)

sparse_matrix_raw <- get_pseudobulked_sc(
  object = sc_object,
  cell_list = pseudobulk_list,
  return_format = "sparse"
)

sparse_matrix_norm <- get_pseudobulked_sc(
  object = sc_object,
  cell_list = pseudobulk_list,
  return_format = "sparse",
  assay = "norm"
)

expect_true(
  current = checkmate::testMatrix(
    matrix_raw,
    mode = "numeric",
    nrows = 3,
    ncols = 81,
    row.names = "named",
    col.names = "named"
  ),
  info = "pseudo-bulking - correct matrix returned: raw and dense"
)

expect_true(
  current = all(round(matrix_raw) == matrix_raw),
  info = "pseudo-bulking: raw counts ARE raw counts"
)

expect_true(
  current = checkmate::testMatrix(
    matrix_raw,
    mode = "numeric",
    nrows = 3,
    ncols = 81,
    row.names = "named",
    col.names = "named"
  ),
  info = "pseudo-bulking - correct matrix returned: raw and norm"
)

expect_true(
  current = checkmate::testClass(sparse_matrix_raw, "dgRMatrix"),
  info = "pseudo-bulking - sparse matrices returned - correct class (1)"
)

expect_true(
  current = checkmate::testClass(sparse_matrix_norm, "dgRMatrix"),
  info = "pseudo-bulking - sparse matrices returned - correct class (2)"
)

expect_equal(
  current = dim(sparse_matrix_raw),
  target = c(3, 81),
  info = "pseudo-bulking - sparse matrices returned - correct dim (1)"
)

expect_equal(
  current = dim(sparse_matrix_norm),
  target = c(3, 81),
  info = "pseudo-bulking - sparse matrices returned - correct dim (2)"
)

expect_true(
  current = all(round(sparse_matrix_raw) == sparse_matrix_raw),
  info = "pseudo-bulking: raw counts ARE raw counts (sparse)"
)

## meta cells ------------------------------------------------------------------

### hdwgcna --------------------------------------------------------------------

hdwgcna <- generate_bt_meta_cells_sc(
  sc_object,
  sc_meta_cell_params = params_sc_bt_metacells(
    target_no_metacells = 50L
  ),
  .verbose = FALSE
)

hdwgcna <- calc_meta_cell_purity(
  hdwgcna,
  original_cell_type = as.character(get_sc_obs(sc_object)$cell_grp)
)

expect_true(
  current = checkmate::testDataTable(
    hdwgcna[[]],
    nrows = 50
  ),
  info = "hdwgcna meta cells obs correct - correct dimensions and type"
)

expect_true(
  current = checkmate::testNames(
    names(hdwgcna[[]]),
    must.include = c(
      "meta_cell_idx",
      "meta_cell_id",
      "no_originating_cells",
      "original_cell_idx"
    )
  ),
  info = "hdwgcna meta cells obs correct - correct columns"
)


expect_true(
  current = mean(hdwgcna[[]]$mc_purity) > 0.75,
  info = "similar cell types are being pulled together"
)

expect_true(
  current = checkmate::testDataTable(
    get_sc_var(hdwgcna)
  ),
  info = "hdwgcna meta cells var correct"
)

expect_equal(
  current = dim(hdwgcna[]),
  target = c(50, 81),
  info = "hdwgcna meta cells - correct return dimensions for the raw counts"
)

expect_true(
  current = checkmate::testClass(hdwgcna[], "dgRMatrix"),
  info = "hdwgcna meta cells - correct compressed sparse matrix type"
)

#### subsetted version ---------------------------------------------------------

# to imitate a scenario in which someone would like to generate meta cells
# in specific cell types; will remove cell_type_3
cells_to_use <- sc_object[[]][cell_grp != "cell_type_3", cell_id]

hdwgcna_small <- generate_bt_meta_cells_sc(
  sc_object,
  sc_meta_cell_params = params_sc_bt_metacells(
    target_no_metacells = 50L,
    max_shared = 10L,
    knn = list(k = 10L)
  ),
  cells_to_use = cells_to_use,
  .verbose = FALSE
)

hdwgcna_small <- calc_meta_cell_purity(
  hdwgcna_small,
  original_cell_type = as.character(get_sc_obs(sc_object)$cell_grp)
)

original_cell_type <- as.character(get_sc_obs(sc_object)$cell_grp)

right_cell_types <- purrr::map_lgl(
  hdwgcna_small[[]]$original_cell_idx,
  function(idx) {
    types <- original_cell_type[idx]
    any(types != "cell_type_3")
  }
)

expect_true(
  current = mean(hdwgcna_small[[]]$mc_purity) > 0.9,
  info = paste(
    "hgwgnca - similar cell types are being pulled together;",
    "subsetted version"
  )
)

expect_true(
  current = all(right_cell_types),
  info = "no unexpected cell types in subsetted version - hdwgcna"
)

### seacells -------------------------------------------------------------------

seacells <- generate_seacells_sc(
  sc_object,
  seacell_params = params_sc_seacells(
    n_sea_cells = 50L,
    min_iter = 5L,
    knn = list(k = 10L)
  ),
  .verbose = FALSE
)

expect_true(
  current = checkmate::testDataTable(
    seacells[[]],
    nrows = 50
  ),
  info = "seacells obs correct - correct dimensions and type"
)

expect_true(
  current = checkmate::testNames(
    names(seacells[[]]),
    must.include = c(
      "meta_cell_idx",
      "meta_cell_id",
      "no_originating_cells",
      "original_cell_idx"
    )
  ),
  info = "seacells obs correct - correct columns"
)

seacells <- calc_meta_cell_purity(
  seacells,
  original_cell_type = as.character(get_sc_obs(sc_object)$cell_grp)
)

expect_true(
  current = mean(seacells[[]]$mc_purity) > 0.9,
  info = "seacell - similar cell types are being pulled together"
)

expect_true(
  current = checkmate::testDataTable(
    get_sc_var(seacells)
  ),
  info = "seacells meta cells var correct"
)

expect_equal(
  current = dim(seacells[]),
  target = c(50, 81),
  info = "seacells - correct return dimensions for the raw counts"
)

expect_true(
  current = checkmate::testClass(seacells[], "dgRMatrix"),
  info = "seacells - correct compressed sparse matrix type"
)

#### smaller subset ------------------------------------------------------------

seacells_small <- generate_seacells_sc(
  sc_object,
  seacell_params = params_sc_seacells(
    n_sea_cells = 50L,
    min_iter = 5L,
    convergence_epsilon = 0.001,
    knn = list(k = 10L)
  ),
  regenerate_knn = TRUE,
  cells_to_use = cells_to_use,
  .verbose = FALSE
)

seacells_small <- calc_meta_cell_purity(
  seacells_small,
  original_cell_type = as.character(get_sc_obs(sc_object)$cell_grp)
)

right_cell_types <- purrr::map_lgl(
  seacells_small[[]]$original_cell_idx,
  function(idx) {
    types <- original_cell_type[idx]
    any(types != "cell_type_3")
  }
)

expect_true(
  current = mean(seacells_small[[]]$mc_purity) > 0.9,
  info = paste(
    "hgwgnca - similar cell types are being pulled together;",
    "subsetted version"
  )
)

expect_true(
  current = all(right_cell_types),
  info = "no unexpected cell types in subsetted version - seacell"
)

### supercells -----------------------------------------------------------------

graining_factor <- 20

supercells <- generate_supercells_sc(
  sc_object,
  sc_supercell_params = params_sc_supercell(
    graining_factor = graining_factor,
    knn = list(k = 10L)
  ),
  regenerate_knn = TRUE,
  .verbose = FALSE
)

expected_meta_cells <- ceiling(sc_object@dims[1] / graining_factor)

expect_true(
  current = checkmate::testDataTable(
    supercells[[]],
    nrows = expected_meta_cells
  ),
  info = "supercells obs correct - correct dimensions and type"
)

expect_true(
  current = checkmate::testNames(
    names(seacells[[]]),
    must.include = c(
      "meta_cell_idx",
      "meta_cell_id",
      "no_originating_cells",
      "original_cell_idx"
    )
  ),
  info = "supercells obs correct - correct columns"
)

supercells <- calc_meta_cell_purity(
  supercells,
  original_cell_type = as.character(get_sc_obs(sc_object)$cell_grp)
)

expect_true(
  current = mean(supercells[[]]$mc_purity) > 0.85,
  info = "supercell - similar cell types are being pulled together"
)

expect_true(
  current = checkmate::testDataTable(
    get_sc_var(supercells)
  ),
  info = "supercell meta cells var correct"
)

expect_equal(
  current = dim(supercells[]),
  target = c(expected_meta_cells, 81),
  info = "supercell - correct return dimensions for the raw counts"
)

expect_true(
  current = checkmate::testClass(supercells[], "dgRMatrix"),
  info = "supercell - correct compressed sparse matrix type"
)

#### smaller subset ------------------------------------------------------------

supercell_small <- generate_supercells_sc(
  sc_object,
  sc_supercell_params = params_sc_supercell(
    graining_factor = graining_factor,
    knn = list(k = 10L)
  ),
  cells_to_use = cells_to_use,
  .verbose = FALSE
)

supercell_small <- calc_meta_cell_purity(
  supercell_small,
  original_cell_type = as.character(get_sc_obs(sc_object)$cell_grp)
)

right_cell_types <- purrr::map_lgl(
  supercell_small[[]]$original_cell_idx,
  function(idx) {
    types <- original_cell_type[idx]
    any(types != "cell_type_3")
  }
)

expect_true(
  current = mean(supercell_small[[]]$mc_purity) > 0.75,
  info = paste(
    "supercell - similar cell types are being pulled together;",
    "subsetted version"
  )
)

expect_true(
  current = all(right_cell_types),
  info = "no unexpected cell types in subsetted version - supercell"
)

## non-contiguous cells to keep -------------------------------------------------

# The tests above load with QC parameters, which drops failing cells at write
# time and leaves cells_to_keep as 0..n. That hides any confusion between the
# embedding row space and the count file row space. Here the data is loaded
# unfiltered and narrowed afterwards via set_cells_to_keep, so cells_to_keep is
# deliberately non-contiguous and the two spaces genuinely differ.

gap_temp_dir <- sc_test_dir("sc_aggregations_gaps")

gap_object <- sc_test_object(
  gap_temp_dir,
  single_cell_test_data,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 1L,
    min_lib_size = 1L,
    min_cells = 1L
  )
)

gap_all_ids <- get_sc_obs(gap_object)$cell_id
gap_n_total <- length(gap_all_ids)

# every other cell, so no working row lines up with its count file row
gap_object <- set_cells_to_keep(
  gap_object,
  gap_all_ids[seq(1L, gap_n_total, by = 2L)]
)

gap_object <- find_hvg_sc(gap_object, hvg_no = hvg_to_keep, .verbose = FALSE)
gap_object <- calculate_pca_sc(gap_object, no_pcs = no_pcs, .verbose = FALSE)
gap_object <- find_neighbours_sc(
  gap_object,
  neighbours_params = params_sc_neighbours(knn = list(k = 15L)),
  .verbose = FALSE
)

gap_keep_1idx <- get_cells_to_keep(gap_object) + 1L

expect_false(
  current = identical(gap_keep_1idx, seq_along(gap_keep_1idx)),
  info = "gap fixture - cells_to_keep is genuinely non-contiguous"
)

# maps a count file row onto the row of the source matrix it was written from
gap_file_to_src <- match(
  get_sc_obs(gap_object)$cell_id,
  rownames(single_cell_test_data$counts)
)

# Aggregated counts must equal the sum of the member cells' raw counts. This is
# the check that actually catches a wrong index space, since misindexed members
# still look plausible on their own.
expect_counts_aggregate <- function(mc_obj, label) {
  members <- mc_obj[[]]$original_cell_idx
  raw <- get_sc_counts(mc_obj, assay = "raw")
  gene_match <- match(colnames(raw), colnames(single_cell_test_data$counts))

  to_check <- seq_len(min(5L, nrow(raw)))
  agree <- purrr::map_lgl(to_check, function(i) {
    src_rows <- gap_file_to_src[members[[i]]]
    manual <- Matrix::colSums(
      single_cell_test_data$counts[src_rows, , drop = FALSE]
    )
    isTRUE(all.equal(as.numeric(manual[gene_match]), as.numeric(raw[i, ])))
  })

  expect_true(
    current = all(unlist(members) %in% gap_keep_1idx),
    info = sprintf("%s - members stay inside cells_to_keep", label)
  )

  expect_true(
    current = all(agree),
    info = sprintf("%s - aggregated counts match the member cells", label)
  )

  expect_equal(
    current = mc_obj@original_assignment$n_cells,
    target = gap_n_total,
    info = sprintf("%s - assignments span the full count file", label)
  )
}

gap_hdwgcna <- generate_bt_meta_cells_sc(
  gap_object,
  sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 50L),
  .verbose = FALSE
)

expect_counts_aggregate(gap_hdwgcna, "hdwgcna non-contiguous")

gap_seacells <- generate_seacells_sc(
  gap_object,
  seacell_params = params_sc_seacells(n_sea_cells = 10L),
  .verbose = FALSE
)

expect_counts_aggregate(gap_seacells, "seacells non-contiguous")

gap_supercells <- generate_supercells_sc(
  gap_object,
  sc_supercell_params = params_sc_supercell(knn = list(k = 10L)),
  .verbose = FALSE
)

expect_counts_aggregate(gap_supercells, "supercell non-contiguous")

### cells_to_use agrees with the default path ----------------------------------

# Asking for every kept cell explicitly must land in the same index space as
# asking for nothing at all.
gap_hdwgcna_use <- generate_bt_meta_cells_sc(
  gap_object,
  sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 50L),
  cells_to_use = get_sc_obs(gap_object, filtered = TRUE)$cell_id,
  .verbose = FALSE
)

expect_true(
  current = all(
    unlist(gap_hdwgcna_use[[]]$original_cell_idx) %in% gap_keep_1idx
  ),
  info = "hdwgcna - cells_to_use path uses the same index space as the default"
)

expect_equal(
  current = gap_hdwgcna_use@original_assignment$n_cells,
  target = gap_n_total,
  info = "hdwgcna - cells_to_use path reports the full count file size"
)

### purity resolves either label vector ----------------------------------------

# The source `cells_to_keep` is recorded on the object, so a QC-passing label
# vector is resolved into the filtered row space instead of mis-attributing cell
# types. Both spaces address the same cells, so the purities must agree.
expect_equal(
  current = calc_meta_cell_purity(
    gap_hdwgcna,
    original_cell_type = unlist(gap_object[["cell_grp"]])
  )[["mc_purity"]]$mc_purity,
  target = calc_meta_cell_purity(
    gap_hdwgcna,
    original_cell_type = as.character(get_sc_obs(gap_object)$cell_grp)
  )[["mc_purity"]]$mc_purity,
  info = "calc_meta_cell_purity agrees on the filtered and unfiltered vectors"
)

expect_silent(
  current = calc_meta_cell_purity(
    gap_hdwgcna,
    original_cell_type = as.character(get_sc_obs(gap_object)$cell_grp)
  ),
  info = "calc_meta_cell_purity accepts the unfiltered label vector"
)

expect_error(
  current = calc_meta_cell_purity(
    gap_hdwgcna,
    original_cell_type = as.character(get_sc_obs(gap_object)$cell_grp)[1:5]
  ),
  info = "calc_meta_cell_purity rejects a vector matching neither index space"
)

### purity optional label and entropy columns ----------------------------------

gap_labels <- as.character(get_sc_obs(gap_object)$cell_grp)

purity_default <- calc_meta_cell_purity(
  gap_hdwgcna,
  original_cell_type = gap_labels
)[[]]

expect_true(
  current = !any(
    c("mc_top_label", "mc_second_label", "mc_second_frac", "mc_entropy") %in%
      colnames(purity_default)
  ),
  info = "calc_meta_cell_purity defaults to the purity column only"
)

purity_top <- calc_meta_cell_purity(
  gap_hdwgcna,
  original_cell_type = gap_labels,
  add_additional_info = "top_label"
)[[]]

expect_true(
  current = all(purity_top$mc_top_label %in% gap_labels),
  info = "calc_meta_cell_purity - the top label is one of the source labels"
)

# the fraction of the reported top label must be exactly the purity
top_label_frac <- purrr::map2_dbl(
  purity_top$original_cell_idx,
  purity_top$mc_top_label,
  \(idx, label) sum(gap_labels[idx] == label) / length(idx)
)

expect_equal(
  current = top_label_frac,
  target = purity_top$mc_purity,
  info = "calc_meta_cell_purity - the top label carries the purity fraction"
)

purity_two <- calc_meta_cell_purity(
  gap_hdwgcna,
  original_cell_type = gap_labels,
  add_additional_info = "top_two_labels",
  add_entropy = TRUE
)[[]]

expect_true(
  current = all(purity_two$mc_second_frac <= purity_two$mc_purity),
  info = "calc_meta_cell_purity - the second label never beats the first"
)

expect_equal(
  current = is.na(purity_two$mc_second_label),
  target = purity_two$mc_purity == 1,
  info = "calc_meta_cell_purity - a pure meta cell has no second label"
)

expect_true(
  current = all(purity_two$mc_entropy >= 0 & purity_two$mc_entropy <= 1),
  info = "calc_meta_cell_purity - the normalised entropy is bounded by [0, 1]"
)

expect_equal(
  current = purity_two$mc_entropy == 0,
  target = purity_two$mc_purity == 1,
  info = "calc_meta_cell_purity - only a pure meta cell has zero entropy"
)

expect_equal(
  current = purity_two$mc_purity,
  target = purity_default$mc_purity,
  info = "calc_meta_cell_purity - the purity is unaffected by the extra columns"
)

expect_error(
  current = calc_meta_cell_purity(
    gap_hdwgcna,
    original_cell_type = gap_labels,
    add_additional_info = "top_three_labels"
  ),
  info = "calc_meta_cell_purity rejects an unknown add_additional_info"
)

## meta cells on a subset -------------------------------------------------------

subset_object <- SingleCellsSubset(
  sc_object = sc_object,
  grouping_column = "cell_grp",
  group = "cell_type_1"
)

subset_object <- find_hvg_sc(
  subset_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)
subset_object <- calculate_pca_sc(
  subset_object,
  no_pcs = no_pcs,
  .verbose = FALSE
)
subset_object <- find_neighbours_sc(
  subset_object,
  neighbours_params = params_sc_neighbours(knn = list(k = 10L)),
  .verbose = FALSE
)

subset_parent_rows <- get_cells_to_keep(subset_object) + 1L
parent_cell_grp <- as.character(get_sc_obs(sc_object)$cell_grp)

expect_subset_metacells <- function(mc_obj, label) {
  members <- unlist(mc_obj[[]]$original_cell_idx)

  expect_true(
    current = all(members %in% subset_parent_rows),
    info = sprintf("%s - members stay inside the subset", label)
  )

  expect_true(
    current = all(parent_cell_grp[members] == "cell_type_1"),
    info = sprintf("%s - members carry the subset group label", label)
  )

  expect_true(
    current = checkmate::testDataTable(get_sc_var(mc_obj)),
    info = sprintf("%s - var table survives the subset", label)
  )
}

subset_hdwgcna <- generate_bt_meta_cells_sc(
  subset_object,
  sc_meta_cell_params = params_sc_bt_metacells(
    target_no_metacells = 20L,
    knn = list(k = 10L)
  ),
  .verbose = FALSE
)

expect_subset_metacells(subset_hdwgcna, "hdwgcna subset")

subset_seacells <- generate_seacells_sc(
  subset_object,
  seacell_params = params_sc_seacells(n_sea_cells = 8L),
  .verbose = FALSE
)

expect_subset_metacells(subset_seacells, "seacells subset")

subset_supercells <- generate_supercells_sc(
  subset_object,
  sc_supercell_params = params_sc_supercell(knn = list(k = 10L)),
  .verbose = FALSE
)

expect_subset_metacells(subset_supercells, "supercell subset")

### pseudo bulking on a subset -------------------------------------------------

subset_cell_ids <- get_sc_obs(subset_object)$cell_id

subset_pseudobulk <- get_pseudobulked_sc(
  object = subset_object,
  cell_list = list(
    grp_a = subset_cell_ids[1:20],
    grp_b = subset_cell_ids[21:40]
  ),
  return_format = "dense",
  .verbose = FALSE
)

expect_equal(
  current = dim(subset_pseudobulk),
  target = c(2L, nrow(get_sc_var(subset_object))),
  info = "pseudo bulking on a subset returns groups x genes"
)

expect_equal(
  current = rownames(subset_pseudobulk),
  target = c("grp_a", "grp_b"),
  info = "pseudo bulking on a subset keeps the group names"
)

# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir, gap_temp_dir)
