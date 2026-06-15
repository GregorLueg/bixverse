# tenx h5 io -------------------------------------------------------------------

library(magrittr)

test_temp_dir <- file.path(tempdir(), "io_tenx")

dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## parameters ------------------------------------------------------------------

min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L

## synthetic data --------------------------------------------------------------

rna <- generate_single_cell_test_data()
adt <- generate_single_cell_test_data_adt()

stopifnot(all(rownames(rna$counts) == rownames(adt$counts)))

### v2 file (single modality) --------------------------------------------------

f_path_v2 <- file.path(test_temp_dir, "tenx_v2.h5")

write_tenx_h5_sc(
  f_path = f_path_v2,
  counts = rna$counts,
  barcodes = rna$obs$cell_id,
  features = data.table::data.table(
    id = rna$var$gene_id,
    name = rna$var$ensembl_id
  ),
  version = "v2"
)

### v3 file (single modality) --------------------------------------------------

f_path_v3 <- file.path(test_temp_dir, "tenx_v3.h5")

write_tenx_h5_sc(
  f_path = f_path_v3,
  counts = rna$counts,
  barcodes = rna$obs$cell_id,
  features = data.table::data.table(
    id = rna$var$gene_id,
    name = rna$var$ensembl_id,
    feature_type = "Gene Expression"
  ),
  version = "v3"
)

### v3 file (multi modal) ------------------------------------------------------

adt_sparse <- as(adt$counts, "RsparseMatrix")
combined_counts <- cbind(rna$counts, adt_sparse)

combined_features <- data.table::data.table(
  id = c(rna$var$gene_id, colnames(adt$counts)),
  name = c(rna$var$ensembl_id, colnames(adt$counts)),
  feature_type = c(
    rep("Gene Expression", nrow(rna$var)),
    rep("Antibody Capture", ncol(adt$counts))
  )
)

f_path_v3_mm <- file.path(test_temp_dir, "tenx_v3_mm.h5")

write_tenx_h5_sc(
  f_path = f_path_v3_mm,
  counts = combined_counts,
  barcodes = rna$obs$cell_id,
  features = combined_features,
  version = "v3"
)

# tests ------------------------------------------------------------------------

## format identification -------------------------------------------------------

meta_v2 <- bixverse:::get_tenx_h5_metadata(f_path_v2)
meta_v3 <- bixverse:::get_tenx_h5_metadata(f_path_v3)
meta_mm <- bixverse:::get_tenx_h5_metadata(f_path_v3_mm)

expect_equal(
  current = meta_v2$version,
  target = "v2",
  info = "tenx h5 - v2 version detection"
)

expect_equal(
  current = meta_v3$version,
  target = "v3",
  info = "tenx h5 - v3 version detection"
)

expect_equal(
  current = c(meta_v2$n_cells, meta_v2$n_genes),
  target = c(nrow(rna$counts), ncol(rna$counts)),
  info = "tenx h5 - v2 dimensions"
)

expect_equal(
  current = c(meta_v3$n_cells, meta_v3$n_genes),
  target = c(nrow(rna$counts), ncol(rna$counts)),
  info = "tenx h5 - v3 dimensions"
)

expect_equal(
  current = meta_mm$n_genes,
  target = ncol(rna$counts) + ncol(adt$counts),
  info = "tenx h5 - multi modal includes all features"
)

## eda function ----------------------------------------------------------------

read_v2 <- read_tenx_h5_metadata(f_path_v2)
read_v3 <- read_tenx_h5_metadata(f_path_v3)
read_mm <- read_tenx_h5_metadata(f_path_v3_mm)

expect_true(
  current = is.null(read_v2$feature_types),
  info = "tenx eda - v2 has no feature_types"
)

expect_equal(
  current = as.integer(read_v3$feature_types["Gene Expression"]),
  target = ncol(rna$counts),
  info = "tenx eda - v3 GEX count correct"
)

expect_equal(
  current = as.integer(read_mm$feature_types["Gene Expression"]),
  target = ncol(rna$counts),
  info = "tenx eda - multi modal GEX count correct"
)

expect_equal(
  current = as.integer(read_mm$feature_types["Antibody Capture"]),
  target = ncol(adt$counts),
  info = "tenx eda - multi modal ADT count correct"
)

expect_equal(
  current = nrow(read_mm$var),
  target = ncol(rna$counts) + ncol(adt$counts),
  info = "tenx eda - var has one row per feature"
)

## expected pass list ----------------------------------------------------------

genes_pass <- which(
  Matrix::colSums(rna$counts != 0) >= min_cells_exp
)

cells_pass <- which(
  (Matrix::rowSums(rna$counts[, genes_pass]) >= min_lib_size) &
    (Matrix::rowSums(rna$counts[, genes_pass] != 0) >= min_genes_exp)
)

expect_true(
  current = length(genes_pass) > 80 & length(genes_pass) != 100,
  info = "tenx - sensible amount of genes pass"
)

expect_true(
  current = length(cells_pass) > 800 & length(cells_pass) != 1000,
  info = "tenx - sensible amount of cells pass"
)

counts_filtered <- rna$counts[cells_pass, genes_pass]
counts_filtered_csc <- as(counts_filtered, "CsparseMatrix")

sc_qc_param <- params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp
)

## direct rust ingestion -------------------------------------------------------

### v2 -------------------------------------------------------------------------

dir_v2 <- file.path(test_temp_dir, "sc_v2")
dir.create(dir_v2, recursive = TRUE, showWarnings = FALSE)

sc_v2 <- SingleCells(dir_data = dir_v2)
rust_v2 <- get_sc_rust_ptr(sc_v2)

file_res_v2 <- rust_v2$tenx_h5_to_file_streaming(
  h5_path = path.expand(f_path_v2),
  version = "v2",
  no_cells = meta_v2$n_cells,
  no_genes = meta_v2$n_genes,
  qc_params = sc_qc_param,
  feature_type = NULL,
  verbose = FALSE
)

expect_equivalent(
  current = file_res_v2$cell_indices + 1,
  target = cells_pass,
  info = "tenx v2 - correct cells being kept"
)

expect_equivalent(
  current = file_res_v2$gene_indices + 1,
  target = genes_pass,
  info = "tenx v2 - correct genes being kept"
)

### v3 -------------------------------------------------------------------------

dir_v3 <- file.path(test_temp_dir, "sc_v3")
dir.create(dir_v3, recursive = TRUE, showWarnings = FALSE)

sc_v3 <- SingleCells(dir_data = dir_v3)
rust_v3 <- get_sc_rust_ptr(sc_v3)

file_res_v3 <- rust_v3$tenx_h5_to_file_streaming(
  h5_path = path.expand(f_path_v3),
  version = "v3",
  no_cells = meta_v3$n_cells,
  no_genes = meta_v3$n_genes,
  qc_params = sc_qc_param,
  feature_type = "Gene Expression",
  verbose = FALSE
)

expect_equivalent(
  current = file_res_v3$cell_indices + 1,
  target = cells_pass,
  info = "tenx v3 - correct cells being kept"
)

expect_equivalent(
  current = file_res_v3$gene_indices + 1,
  target = genes_pass,
  info = "tenx v3 - correct genes being kept"
)

### auto version detection -----------------------------------------------------

dir_auto <- file.path(test_temp_dir, "sc_auto")
dir.create(dir_auto, recursive = TRUE, showWarnings = FALSE)

sc_auto <- SingleCells(dir_data = dir_auto)
rust_auto <- get_sc_rust_ptr(sc_auto)

file_res_auto <- rust_auto$tenx_h5_to_file_streaming(
  h5_path = path.expand(f_path_v3),
  version = "auto",
  no_cells = meta_v3$n_cells,
  no_genes = meta_v3$n_genes,
  qc_params = sc_qc_param,
  feature_type = "Gene Expression",
  verbose = FALSE
)

expect_equivalent(
  current = file_res_auto$cell_indices + 1,
  target = cells_pass,
  info = "tenx - auto version detection - cells"
)

expect_equivalent(
  current = file_res_auto$gene_indices + 1,
  target = genes_pass,
  info = "tenx - auto version detection - genes"
)

### multi modal: only GEX kept -------------------------------------------------

dir_mm <- file.path(test_temp_dir, "sc_mm_rust")
dir.create(dir_mm, recursive = TRUE, showWarnings = FALSE)

sc_mm_rust <- SingleCells(dir_data = dir_mm)
rust_mm <- get_sc_rust_ptr(sc_mm_rust)

file_res_mm <- rust_mm$tenx_h5_to_file_streaming(
  h5_path = path.expand(f_path_v3_mm),
  version = "v3",
  no_cells = meta_mm$n_cells,
  no_genes = meta_mm$n_genes,
  qc_params = sc_qc_param,
  feature_type = "Gene Expression",
  verbose = FALSE
)

expect_equivalent(
  current = file_res_mm$cell_indices + 1,
  target = cells_pass,
  info = "tenx multi modal - same cells as RNA only"
)

expect_equivalent(
  current = file_res_mm$gene_indices + 1,
  target = genes_pass,
  info = "tenx multi modal - only GEX features kept"
)

expect_true(
  current = max(file_res_mm$gene_indices + 1) <= ncol(rna$counts),
  info = "tenx multi modal - no ADT index leaked into gene_indices"
)

### edge case: absent feature_type ---------------------------------------------

expect_error(
  current = rust_mm$tenx_h5_to_file_streaming(
    h5_path = path.expand(f_path_v3_mm),
    version = "v3",
    no_cells = meta_mm$n_cells,
    no_genes = meta_mm$n_genes,
    qc_params = sc_qc_param,
    feature_type = "Surface Protein",
    verbose = FALSE
  ),
  info = "tenx multi modal - requesting absent feature_type errors"
)

## full object load ------------------------------------------------------------

dir_full <- file.path(test_temp_dir, "sc_full")
dir.create(dir_full, recursive = TRUE, showWarnings = FALSE)

sc_full <- SingleCells(dir_data = dir_full)

sc_full <- load_tenx_h5(
  object = sc_full,
  h5_path = f_path_v3,
  sc_qc_param = sc_qc_param,
  feature_type = "Gene Expression",
  streaming = 0L,
  .verbose = FALSE
)

### obs ------------------------------------------------------------------------

obs_object <- sc_full[[]]

expect_equal(
  current = obs_object$cell_id,
  target = rna$obs$cell_id[cells_pass],
  info = "tenx v3 (full load) - correct cells kept"
)

expect_equivalent(
  current = Matrix::rowSums(counts_filtered),
  target = obs_object$lib_size,
  info = "tenx v3 (full load) - lib size correct"
)

expect_equivalent(
  current = Matrix::rowSums(counts_filtered != 0),
  target = obs_object$nnz,
  info = "tenx v3 (full load) - nnz correct"
)

### var ------------------------------------------------------------------------

vars_object <- get_sc_var(sc_full)

expect_equal(
  current = vars_object$gene_id,
  target = rna$var$gene_id[genes_pass],
  info = "tenx v3 (full load) - correct genes kept"
)

expect_true(
  current = "gene_name" %in% names(vars_object),
  info = "tenx v3 (full load) - gene_name populated"
)

expect_true(
  current = "feature_type" %in% names(vars_object),
  info = "tenx v3 (full load) - feature_type populated"
)

expect_true(
  current = all(vars_object$feature_type == "Gene Expression"),
  info = "tenx v3 (full load) - only GEX surfaced"
)

expect_true(
  current = checkmate::qtest(vars_object$no_cells_exp, "I+"),
  info = "tenx v3 (full load) - gene NNZ populated"
)

### counts ---------------------------------------------------------------------

counts_csr_obj <- sc_full[,, return_format = "cell"]

expect_equal(
  current = counts_csr_obj,
  target = counts_filtered,
  info = "tenx v3 (full load) - full CSR retrieval"
)

set.seed(42L)
random_cells <- sample(1:nrow(counts_csr_obj), size = 20L)

expect_equal(
  current = sc_full[random_cells, , return_format = "cell"],
  target = counts_filtered[random_cells, ],
  info = "tenx v3 (full load) - CSR random subset"
)

counts_csc_obj <- sc_full[,, return_format = "gene"]

expect_equal(
  current = counts_csc_obj,
  target = counts_filtered_csc,
  info = "tenx v3 (full load) - full CSC retrieval"
)

set.seed(123L)
random_genes <- sample(1:ncol(counts_csr_obj), size = 20L)

expect_equal(
  current = sc_full[, random_genes, return_format = "gene"],
  target = counts_filtered_csc[, random_genes],
  info = "tenx v3 (full load) - CSC random subset"
)

### getters --------------------------------------------------------------------

expect_equal(
  current = get_cell_names(sc_full),
  target = rna$obs$cell_id[cells_pass],
  info = "tenx v3 (full load) - cell names"
)

expect_equal(
  current = get_gene_names(sc_full),
  target = rna$var$gene_id[genes_pass],
  info = "tenx v3 (full load) - gene names"
)

## v2 full object load ---------------------------------------------------------

dir_full_v2 <- file.path(test_temp_dir, "sc_full_v2")
dir.create(dir_full_v2, recursive = TRUE, showWarnings = FALSE)

sc_full_v2 <- SingleCells(dir_data = dir_full_v2)

sc_full_v2 <- load_tenx_h5(
  object = sc_full_v2,
  h5_path = f_path_v2,
  sc_qc_param = sc_qc_param,
  streaming = 0L,
  .verbose = FALSE
)

expect_equal(
  current = sc_full_v2[[]]$cell_id,
  target = rna$obs$cell_id[cells_pass],
  info = "tenx v2 (full load) - correct cells"
)

expect_equal(
  current = get_sc_var(sc_full_v2)$gene_id,
  target = rna$var$gene_id[genes_pass],
  info = "tenx v2 (full load) - correct genes"
)

expect_equal(
  current = sc_full_v2[,, return_format = "cell"],
  target = counts_filtered,
  info = "tenx v2 (full load) - counts"
)

## ADT reading via read_tenx_h5_adt --------------------------------------------

adt_read <- read_tenx_h5_adt(
  f_path = f_path_v3_mm,
  feature_type = "Antibody Capture"
)

expect_true(
  current = checkmate::testMatrix(
    adt_read,
    row.names = "named",
    col.names = "named"
  ),
  info = "tenx ADT read - returns a named matrix"
)

expect_equal(
  current = dim(adt_read),
  target = c(nrow(rna$counts), ncol(adt$counts)),
  info = "tenx ADT read - correct dimensions"
)

expect_equal(
  current = rownames(adt_read),
  target = rna$obs$cell_id,
  info = "tenx ADT read - correct barcodes"
)

expect_equal(
  current = colnames(adt_read),
  target = colnames(adt$counts),
  info = "tenx ADT read - correct feature names"
)

expect_equivalent(
  current = unname(adt_read),
  target = unname(adt$counts),
  info = "tenx ADT read - values match input"
)

## multi modal end-to-end ------------------------------------------------------

dir_mm_full <- file.path(test_temp_dir, "sc_mm_full")
dir.create(dir_mm_full, recursive = TRUE, showWarnings = FALSE)

sc_mm <- SingleCellsMultiModal(dir_data = dir_mm_full)

sc_mm <- load_tenx_h5(
  object = sc_mm,
  h5_path = f_path_v3_mm,
  sc_qc_param = sc_qc_param,
  feature_type = "Gene Expression",
  streaming = 0L,
  .verbose = FALSE
)

expect_equal(
  current = length(get_cells_to_keep(sc_mm)),
  target = length(cells_pass),
  info = "tenx multi modal (full) - same cells as RNA only"
)

expect_equal(
  current = get_gene_names(sc_mm),
  target = rna$var$gene_id[genes_pass],
  info = "tenx multi modal (full) - only GEX genes in var table"
)

adt_for_class <- read_tenx_h5_adt(
  f_path = f_path_v3_mm,
  feature_type = "Antibody Capture"
)

sc_mm <- add_adt_counts_sc(
  sc_mm,
  adt_counts = adt_for_class,
  method = "clr"
)

expect_true(
  current = !is.null(S7::prop(sc_mm, "adt_counts")),
  info = "tenx multi modal - ADT attached"
)

expect_equal(
  current = get_adt_names(sc_mm),
  target = colnames(adt$counts),
  info = "tenx multi modal - ADT feature names preserved"
)

adt_in_obj <- get_sc_counts(
  sc_mm,
  modality = "adt",
  assay = "raw",
  .verbose = FALSE
)

expect_equal(
  current = nrow(adt_in_obj),
  target = length(cells_pass),
  info = "tenx multi modal - ADT subset to kept cells"
)

expect_equivalent(
  current = unname(adt_in_obj),
  target = unname(adt$counts[cells_pass, ]),
  info = "tenx multi modal - ADT values match input on kept cells"
)

## multi file ingestion -------------------------------------------------------

min_lib_size <- 300L
min_genes_exp <- 25L
min_cells_exp <- 250L

sc_qc_param <- params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp
)

### synthetic data ------------------------------------------------------------

# split the existing rna data into two disjoint cell sets with overlapping
# but non-identical gene sets to exercise the universe logic
set.seed(7L)
cell_perm <- sample(seq_len(nrow(rna$counts)))
half_n <- length(cell_perm) %/% 2L
cells_a <- sort(cell_perm[1:half_n])
cells_b <- sort(cell_perm[(half_n + 1L):length(cell_perm)])

genes_a <- 1:80
genes_b <- 21:100
genes_isect <- intersect(genes_a, genes_b)
genes_union_idx <- union(genes_a, genes_b)

counts_a <- rna$counts[cells_a, genes_a]
counts_b <- rna$counts[cells_b, genes_b]

barcodes_a <- paste0("A_", rna$obs$cell_id[cells_a])
barcodes_b <- paste0("B_", rna$obs$cell_id[cells_b])

features_a_v3 <- data.table::data.table(
  id = rna$var$gene_id[genes_a],
  name = rna$var$ensembl_id[genes_a],
  feature_type = "Gene Expression"
)
features_b_v3 <- data.table::data.table(
  id = rna$var$gene_id[genes_b],
  name = rna$var$ensembl_id[genes_b],
  feature_type = "Gene Expression"
)

f_path_multi_a <- file.path(test_temp_dir, "multi_a.h5")
f_path_multi_b <- file.path(test_temp_dir, "multi_b.h5")

write_tenx_h5_sc(
  f_path = f_path_multi_a,
  counts = counts_a,
  barcodes = barcodes_a,
  features = features_a_v3,
  version = "v3"
)

write_tenx_h5_sc(
  f_path = f_path_multi_b,
  counts = counts_b,
  barcodes = barcodes_b,
  features = features_b_v3,
  version = "v3"
)

# mixed version: write file_b as v2 instead
features_b_v2 <- data.table::data.table(
  id = rna$var$gene_id[genes_b],
  name = rna$var$ensembl_id[genes_b]
)
f_path_multi_b_v2 <- file.path(test_temp_dir, "multi_b_v2.h5")
write_tenx_h5_sc(
  f_path = f_path_multi_b_v2,
  counts = counts_b,
  barcodes = barcodes_b,
  features = features_b_v2,
  version = "v2"
)

# multi modal file_b: add some ADT features that must be dropped
adt_sub <- adt$counts[cells_b, ]
adt_sub_sparse <- as(adt_sub, "RsparseMatrix")
counts_b_mm <- cbind(counts_b, adt_sub_sparse)
features_b_mm <- data.table::data.table(
  id = c(rna$var$gene_id[genes_b], colnames(adt$counts)),
  name = c(rna$var$ensembl_id[genes_b], colnames(adt$counts)),
  feature_type = c(
    rep("Gene Expression", length(genes_b)),
    rep("Antibody Capture", ncol(adt$counts))
  )
)
f_path_multi_b_mm <- file.path(test_temp_dir, "multi_b_mm.h5")
write_tenx_h5_sc(
  f_path = f_path_multi_b_mm,
  counts = counts_b_mm,
  barcodes = barcodes_b,
  features = features_b_mm,
  version = "v3"
)

### prescan: intersection ----------------------------------------------------

prescan_isect <- prescan_tenx_h5_files(
  h5_paths = c(A = f_path_multi_a, B = f_path_multi_b),
  feature_type = "Gene Expression",
  gene_universe = "intersection",
  .verbose = FALSE
)

expect_equal(
  current = prescan_isect$universe_size,
  target = length(genes_isect),
  info = "multi tenx prescan - intersection universe size"
)

expect_true(
  current = setequal(
    prescan_isect$universe,
    rna$var$gene_id[genes_isect]
  ),
  info = "multi tenx prescan - intersection universe content"
)

expect_equal(
  current = names(prescan_isect$file_tasks),
  target = c("A", "B"),
  info = "multi tenx prescan - file_tasks names match exp_ids"
)

expect_equal(
  current = prescan_isect$file_tasks[["A"]]$version,
  target = "v3",
  info = "multi tenx prescan - version captured"
)

### prescan: union -----------------------------------------------------------

prescan_union <- prescan_tenx_h5_files(
  h5_paths = c(A = f_path_multi_a, B = f_path_multi_b),
  feature_type = "Gene Expression",
  gene_universe = "union",
  .verbose = FALSE
)

expect_equal(
  current = prescan_union$universe_size,
  target = length(genes_union_idx),
  info = "multi tenx prescan - union universe size"
)

### prescan: ADT filtered out from multi modal file --------------------------

prescan_mm <- prescan_tenx_h5_files(
  h5_paths = c(A = f_path_multi_a, B = f_path_multi_b_mm),
  feature_type = "Gene Expression",
  gene_universe = "intersection",
  .verbose = FALSE
)

expect_equal(
  current = prescan_mm$universe_size,
  target = length(genes_isect),
  info = "multi tenx prescan - ADT features excluded from universe"
)

expect_true(
  current = all(
    is.na(prescan_mm$file_tasks[["B"]]$gene_local_to_universe[
      (length(genes_b) + 1L):(length(genes_b) + ncol(adt$counts))
    ])
  ),
  info = "multi tenx prescan - ADT features map to NA"
)

### prescan: missing feature_type errors -------------------------------------

expect_error(
  current = prescan_tenx_h5_files(
    h5_paths = c(A = f_path_multi_a, B = f_path_multi_b),
    feature_type = "Surface Protein",
    gene_universe = "intersection",
    .verbose = FALSE
  ),
  info = "multi tenx prescan - unknown feature_type errors"
)

### expected pass list (intersection) ----------------------------------------

# global gene QC is across both files
counts_a_isect <- rna$counts[cells_a, genes_isect]
counts_b_isect <- rna$counts[cells_b, genes_isect]
nnz_a <- Matrix::colSums(counts_a_isect != 0)
nnz_b <- Matrix::colSums(counts_b_isect != 0)
genes_pass_global <- which((nnz_a + nnz_b) >= min_cells_exp)

cells_a_pass <- which(
  (Matrix::rowSums(counts_a_isect[, genes_pass_global]) >= min_lib_size) &
    (Matrix::rowSums(counts_a_isect[, genes_pass_global] != 0) >= min_genes_exp)
)
cells_b_pass <- which(
  (Matrix::rowSums(counts_b_isect[, genes_pass_global]) >= min_lib_size) &
    (Matrix::rowSums(counts_b_isect[, genes_pass_global] != 0) >= min_genes_exp)
)

### full multi load: two v3 files --------------------------------------------

dir_multi <- file.path(test_temp_dir, "sc_multi_tenx")
dir.create(dir_multi, recursive = TRUE, showWarnings = FALSE)

sc_multi <- SingleCells(dir_data = dir_multi)

sc_multi <- load_multi_tenx_h5(
  object = sc_multi,
  prescan_result = prescan_isect,
  sc_qc_param = sc_qc_param,
  streaming = 0L,
  .verbose = FALSE
)

expect_equal(
  current = dim(sc_multi)[1],
  target = length(cells_a_pass) + length(cells_b_pass),
  info = "multi tenx (full) - total cells = sum of per-file passes"
)

expect_equal(
  current = dim(sc_multi)[2],
  target = length(genes_pass_global),
  info = "multi tenx (full) - genes match global QC"
)

obs_multi <- sc_multi[[]]

expect_true(
  current = "exp_id" %in% names(obs_multi),
  info = "multi tenx (full) - exp_id column added"
)

expect_equal(
  current = sort(unique(obs_multi$exp_id)),
  target = c("A", "B"),
  info = "multi tenx (full) - both files represented"
)

expect_equal(
  current = sum(obs_multi$exp_id == "A"),
  target = length(cells_a_pass),
  info = "multi tenx (full) - file A cell count"
)

expect_equal(
  current = sum(obs_multi$exp_id == "B"),
  target = length(cells_b_pass),
  info = "multi tenx (full) - file B cell count"
)

expect_true(
  current = all(grepl("^A_", obs_multi[exp_id == "A", cell_id])) &&
    all(grepl("^B_", obs_multi[exp_id == "B", cell_id])),
  info = "multi tenx (full) - cell_id prefixed with exp_id"
)

var_multi <- get_sc_var(sc_multi)
expect_true(
  current = setequal(
    var_multi$gene_id,
    rna$var$gene_id[genes_isect][genes_pass_global]
  ),
  info = "multi tenx (full) - gene ids match QC pass"
)

### count round-trip per file ------------------------------------------------

# rebuild the expected dense counts from the (cells_a x universe-survivors)
# and (cells_b x universe-survivors) slices, in the universe gene order

gene_id_kept <- rna$var$gene_id[genes_isect][genes_pass_global]
universe_order <- match(prescan_isect$universe, rna$var$gene_id)
gene_id_universe_order <- rna$var$gene_id[universe_order]
gene_idx_final_in_universe <- match(gene_id_kept, gene_id_universe_order)

counts_a_expected <- counts_a_isect[cells_a_pass, , drop = FALSE]
colnames(counts_a_expected) <- rna$var$gene_id[genes_isect]
counts_a_expected <- counts_a_expected[, gene_id_kept, drop = FALSE]

counts_b_expected <- counts_b_isect[cells_b_pass, , drop = FALSE]
colnames(counts_b_expected) <- rna$var$gene_id[genes_isect]
counts_b_expected <- counts_b_expected[, gene_id_kept, drop = FALSE]

counts_multi_obj <- sc_multi[,, return_format = "cell"]

expect_equivalent(
  current = unname(as.matrix(
    counts_multi_obj[obs_multi$exp_id == "A", ]
  )),
  target = unname(as.matrix(counts_a_expected)),
  info = "multi tenx (full) - file A counts round-trip"
)

expect_equivalent(
  current = unname(as.matrix(
    counts_multi_obj[obs_multi$exp_id == "B", ]
  )),
  target = unname(as.matrix(counts_b_expected)),
  info = "multi tenx (full) - file B counts round-trip"
)

### mixed v2 + v3 ------------------------------------------------------------

prescan_mixed <- prescan_tenx_h5_files(
  h5_paths = c(A = f_path_multi_a, B = f_path_multi_b_v2),
  feature_type = "Gene Expression",
  gene_universe = "intersection",
  .verbose = FALSE
)

expect_equal(
  current = c(
    prescan_mixed$file_tasks[["A"]]$version,
    prescan_mixed$file_tasks[["B"]]$version
  ),
  target = c("v3", "v2"),
  info = "multi tenx prescan - mixed versions captured"
)

dir_mixed <- file.path(test_temp_dir, "sc_multi_tenx_mixed")
dir.create(dir_mixed, recursive = TRUE, showWarnings = FALSE)
sc_mixed <- SingleCells(dir_data = dir_mixed)

sc_mixed <- load_multi_tenx_h5(
  object = sc_mixed,
  prescan_result = prescan_mixed,
  sc_qc_param = sc_qc_param,
  streaming = 0L,
  .verbose = FALSE
)

expect_equal(
  current = dim(sc_mixed)[1],
  target = length(cells_a_pass) + length(cells_b_pass),
  info = "multi tenx mixed v2/v3 - same total cells as both-v3 case"
)

expect_equal(
  current = dim(sc_mixed)[2],
  target = length(genes_pass_global),
  info = "multi tenx mixed v2/v3 - same genes as both-v3 case"
)

### multi modal: ADT does not leak into counts -------------------------------

dir_mm_multi <- file.path(test_temp_dir, "sc_multi_tenx_mm")
dir.create(dir_mm_multi, recursive = TRUE, showWarnings = FALSE)
sc_mm_multi <- SingleCells(dir_data = dir_mm_multi)

sc_mm_multi <- load_multi_tenx_h5(
  object = sc_mm_multi,
  prescan_result = prescan_mm,
  sc_qc_param = sc_qc_param,
  streaming = 0L,
  .verbose = FALSE
)

expect_equal(
  current = dim(sc_mm_multi)[2],
  target = length(genes_pass_global),
  info = "multi tenx multi-modal - ADT excluded from genes"
)

expect_true(
  current = !any(get_gene_names(sc_mm_multi) %in% colnames(adt$counts)),
  info = "multi tenx multi-modal - no ADT names in gene_id"
)

# clean up ---------------------------------------------------------------------

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
