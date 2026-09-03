# SingleCellExperiment i/o tests -----------------------------------------------

source("helper_sc.R", local = TRUE)

# the Bioconductor classes are Suggests, so this file is a no-op without them
if (
  !requireNamespace("SingleCellExperiment", quietly = TRUE) ||
    !requireNamespace("SummarizedExperiment", quietly = TRUE)
) {
  exit_file("SingleCellExperiment / SummarizedExperiment not installed")
}

test_temp_dir <- sc_test_dir("sc_io_sce")

## synthetic test data ---------------------------------------------------------

set.seed(42L)

n_genes <- 120L
n_cells <- 400L

# dense enough that every gene and cell clears the QC floors below
counts_gc <- Matrix::rsparsematrix(
  nrow = n_genes,
  ncol = n_cells,
  density = 0.4,
  rand.x = \(n) as.numeric(rpois(n, 5) + 1L)
)
counts_gc <- as(counts_gc, "dgCMatrix")

gene_ids <- sprintf("gene_%03i", seq_len(n_genes))
cell_ids <- sprintf("cell_%03i", seq_len(n_cells))
rownames(counts_gc) <- gene_ids
colnames(counts_gc) <- cell_ids

col_data <- data.frame(
  donor = rep(sprintf("donor_%i", 1:8), length.out = n_cells),
  condition = rep(c("ctrl", "stim"), length.out = n_cells),
  stringsAsFactors = FALSE
)
rownames(col_data) <- cell_ids

row_data <- data.frame(
  ensembl = sprintf("ENSG%08i", seq_len(n_genes)),
  chromosome = rep(c("1", "2", "X"), length.out = n_genes),
  stringsAsFactors = FALSE
)
rownames(row_data) <- gene_ids

make_sce <- function(counts = counts_gc, cd = col_data, rd = row_data) {
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = cd,
    rowData = rd
  )
}

# SummarizedExperiment takes the dimnames off the rowData when the matrix has
# none, so both sides have to be stripped to get an object with no gene names
make_sce_unnamed_genes <- function(rd = row_data) {
  counts <- counts_gc
  rownames(counts) <- NULL
  rownames(rd) <- NULL
  make_sce(counts = counts, rd = rd)
}

# permissive, so nothing is filtered and the round trip is exact
qc_open <- params_sc_min_quality(
  min_unique_genes = 1L,
  min_lib_size = 1L,
  min_cells = 1L,
  target_size = 1e4
)

load_it <- function(sce, name, ...) {
  d <- file.path(test_temp_dir, name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  obj <- SingleCells(dir_data = d)
  load_sce(obj, sce = sce, sc_qc_param = qc_open, .verbose = FALSE, ...)
}

# tests ------------------------------------------------------------------------

## the happy path --------------------------------------------------------------

sce <- make_sce()
obj <- load_it(sce, "basic")

expect_equal(
  current = S7::prop(obj, "dims"),
  target = c(n_cells, n_genes),
  info = "load_sce transposes genes x cells into cells x genes"
)

obs <- obj[[]]
var <- get_sc_var(obj)

expect_equal(
  current = obs$cell_id,
  target = cell_ids,
  info = "colnames of the SCE become cell_id, in order"
)

expect_equal(
  current = var$gene_id,
  target = gene_ids,
  info = "rownames of the SCE become gene_id, in order"
)

expect_true(
  current = all(c("donor", "condition") %in% names(obs)),
  info = "colData columns are carried into obs"
)

expect_true(
  current = all(c("ensembl", "chromosome") %in% names(var)),
  info = "rowData columns are carried into var"
)

expect_equal(
  current = obs$donor,
  target = col_data$donor,
  info = "the colData values survive the trip"
)

## the counts survive ----------------------------------------------------------

expect_equal(
  current = obs$lib_size,
  target = as.numeric(Matrix::colSums(counts_gc)),
  info = "library sizes match the source matrix"
)

expect_equal(
  current = obs$nnz,
  target = as.integer(Matrix::colSums(counts_gc > 0)),
  info = "non-zero counts per cell match the source matrix"
)

expect_equal(
  current = var$no_cells_exp,
  target = as.integer(Matrix::rowSums(counts_gc > 0)),
  info = "cells expressing each gene match the source matrix"
)

## it agrees with load_r_data --------------------------------------------------

# load_sce is a shim over load_r_data, so the two have to land in the same place
d_direct <- file.path(test_temp_dir, "direct")
dir.create(d_direct, recursive = TRUE, showWarnings = FALSE)
obj_direct <- load_r_data(
  object = SingleCells(dir_data = d_direct),
  counts = new(
    "dgRMatrix",
    p = counts_gc@p,
    j = counts_gc@i,
    x = counts_gc@x,
    Dim = rev(counts_gc@Dim)
  ),
  obs = data.table::data.table(barcode = cell_ids, col_data),
  var = data.table::data.table(gene_id = gene_ids, row_data),
  sc_qc_param = qc_open,
  .verbose = FALSE
)

expect_equal(
  current = get_sc_var(obj_direct)$gene_id,
  target = var$gene_id,
  info = "load_sce and load_r_data agree on the genes"
)

expect_equal(
  current = obj_direct[[]]$lib_size,
  target = obs$lib_size,
  info = "load_sce and load_r_data agree on the counts"
)

## identifiers absent from the dimnames ----------------------------------------

# Bioconductor objects routinely leave the dimnames empty and keep the
# identifiers in the metadata, e.g. the ageing thymus data
sce_no_rn <- make_sce_unnamed_genes()

expect_warning(
  current = load_it(sce_no_rn, "no_rownames"),
  pattern = "No gene names",
  info = "falling back to the metadata for gene ids is announced"
)

obj_no_rn <- suppressWarnings(load_it(sce_no_rn, "no_rownames2"))

expect_equal(
  current = get_sc_var(obj_no_rn)$gene_id,
  target = row_data$ensembl,
  info = "the first rowData column stands in for absent rownames"
)

## identifiers that are not a usable key ---------------------------------------

row_data_na <- row_data
row_data_na$ensembl[c(3L, 9L, 20L)] <- NA_character_
sce_na <- make_sce_unnamed_genes(rd = row_data_na)

expect_warning(
  current = load_it(sce_na, "na_ids"),
  pattern = "no identifier",
  info = "missing identifiers are announced rather than passed through"
)

obj_na <- suppressWarnings(load_it(sce_na, "na_ids2"))
var_na <- get_sc_var(obj_na)

expect_true(
  current = !anyNA(var_na$gene_id) && !anyDuplicated(var_na$gene_id),
  info = "the repaired gene ids are a usable join key"
)

row_data_dup <- row_data
row_data_dup$ensembl[5L] <- row_data_dup$ensembl[4L]
sce_dup <- make_sce_unnamed_genes(rd = row_data_dup)

expect_warning(
  current = load_it(sce_dup, "dup_ids"),
  pattern = "duplicated",
  info = "duplicate identifiers are announced"
)

obj_dup <- suppressWarnings(load_it(sce_dup, "dup_ids2"))

expect_true(
  current = !anyDuplicated(get_sc_var(obj_dup)$gene_id),
  info = "duplicate gene ids are made unique"
)

## metadata that would collide with the reserved names -------------------------

# the first obs column is forced to cell_id downstream, so a CellID column in
# the colData would otherwise come back as cell_id.1
col_data_clash <- col_data
col_data_clash$CellID <- cell_ids
sce_clash <- make_sce(cd = col_data_clash)

obs_clash <- load_it(sce_clash, "clash")[[]]

expect_false(
  current = "cell_id.1" %in% names(obs_clash),
  info = "a colliding colData column does not produce a suffixed duplicate"
)

expect_equal(
  current = obs_clash$cell_id,
  target = cell_ids,
  info = "the collision is resolved in favour of the real identifiers"
)

## non-atomic metadata ---------------------------------------------------------

sce_nested <- make_sce()
SummarizedExperiment::colData(sce_nested)$nested <- S4Vectors::DataFrame(
  a = seq_len(n_cells)
)

expect_warning(
  current = load_it(sce_nested, "nested"),
  pattern = "non-atomic",
  info = "columns the DuckDB cannot take are dropped with a warning"
)

## input validation ------------------------------------------------------------

expect_error(
  current = load_it(sce, "bad_assay", assay_name = "logcounts"),
  pattern = "not found",
  info = "load_sce rejects an assay the object does not have"
)

sce_neg <- make_sce()
SummarizedExperiment::assay(sce_neg, "counts")[1L, 1L] <- -5

expect_error(
  current = load_it(sce_neg, "negative"),
  pattern = "negative values",
  info = "load_sce rejects an assay that is not raw counts"
)

sce_dense <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = as.matrix(counts_gc)),
  colData = col_data,
  rowData = row_data
)

expect_error(
  current = load_it(sce_dense, "dense"),
  pattern = "not a dgCMatrix",
  info = "load_sce rejects an assay it cannot reinterpret"
)

# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
