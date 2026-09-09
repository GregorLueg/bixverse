# cellsweep --------------------------------------------------------------------

source("helper_sc.R", local = TRUE)

set.seed(123L)

test_dir_input <- sc_test_dir("cellsweep_input")
test_dir_output <- sc_test_dir("cellsweep_output")
test_dir_inferred <- sc_test_dir("cellsweep_inferred")
test_dir_filtered <- sc_test_dir("cellsweep_filtered")

## synthetic data --------------------------------------------------------------

fixture <- cellsweep_test_fixture()

lib_size <- as.integer(Matrix::rowSums(fixture$counts))
is_empty_true <- fixture$obs$is_empty
n_real <- fixture$n_real
n_genes <- ncol(fixture$counts)

# the soup is cell type 1 plus background, so genes 1:20 are cell type 1's
# marker block and anything cell type 2 carries there arrived via the ambient
soup_block <- seq_len(20L)

## helper functions ------------------------------------------------------------

empty_droplet_recall <- function(call) {
  sum(call & is_empty_true) / sum(is_empty_true)
}

empty_droplet_precision <- function(call) {
  sum(call & is_empty_true) / max(sum(call), 1L)
}

soup_fraction <- function(counts) {
  sum(counts[, soup_block]) / sum(counts)
}

# tests ------------------------------------------------------------------------

## empty droplet inference -----------------------------------------------------

### rust logic -----------------------------------------------------------------

umi_call <- bixverse:::rs_sc_infer_empty_droplets(
  lib_size = lib_size,
  empty_params = unclass(
    params_sc_empty_droplets(method = "umi_cutoff", umi_cutoff = 500L)
  )
)

expect_true(
  current = checkmate::qtest(umi_call, sprintf("B%s", length(lib_size))),
  info = "rust empty droplets: umi cutoff returns a logical per barcode"
)

expect_true(
  current = empty_droplet_recall(umi_call) == 1,
  info = "rust empty droplets: umi cutoff finds every empty droplet"
)

expect_true(
  current = empty_droplet_precision(umi_call) == 1,
  info = "rust empty droplets: umi cutoff calls no real barcode empty"
)

expected_call <- bixverse:::rs_sc_infer_empty_droplets(
  lib_size = lib_size,
  empty_params = unclass(
    params_sc_empty_droplets(
      method = "expected_cells",
      expected_cells = n_real
    )
  )
)

expect_true(
  current = empty_droplet_recall(expected_call) == 1,
  info = "rust empty droplets: expected cells finds every empty droplet"
)

expect_true(
  current = sum(!expected_call) == n_real,
  info = "rust empty droplets: expected cells keeps exactly that many barcodes"
)

knee_call <- bixverse:::rs_sc_infer_empty_droplets(
  lib_size = lib_size,
  empty_params = unclass(params_sc_empty_droplets(method = "knee"))
)

expect_true(
  current = empty_droplet_recall(knee_call) == 1,
  info = "rust empty droplets: knee finds every empty droplet"
)

# the knee sits a few ranks ahead of the cliff, so it sweeps up real barcodes
expect_true(
  current = empty_droplet_precision(knee_call) >= 0.99,
  info = "rust empty droplets: knee oversweeps only marginally"
)

### resolve helper -------------------------------------------------------------

obs_resolve <- data.table::data.table(
  lib_size = lib_size,
  is_empty = is_empty_true
)

expect_equal(
  current = bixverse:::.resolve_empty_droplets(
    obs = obs_resolve,
    empty_params = params_sc_empty_droplets(
      method = "supplied",
      is_empty_column = "is_empty"
    ),
    .verbose = FALSE
  ),
  target = is_empty_true,
  info = "resolve empty droplets: supplied hands back the obs column"
)

expect_equal(
  current = bixverse:::.resolve_empty_droplets(
    obs = obs_resolve,
    empty_params = params_sc_empty_droplets(
      method = "umi_cutoff",
      umi_cutoff = 500L
    ),
    .verbose = FALSE
  ),
  target = umi_call,
  info = "resolve empty droplets: inference matches the rust call"
)

### error handling -------------------------------------------------------------

expect_error(
  current = bixverse:::.resolve_empty_droplets(
    obs = obs_resolve,
    empty_params = params_sc_empty_droplets(
      method = "supplied",
      is_empty_column = "not_a_column"
    ),
    .verbose = FALSE
  ),
  pattern = "is not a column in the obs table",
  info = "resolve empty droplets: missing column errors"
)

expect_error(
  current = bixverse:::.resolve_empty_droplets(
    obs = data.table::data.table(
      lib_size = lib_size,
      is_empty = as.integer(is_empty_true)
    ),
    empty_params = params_sc_empty_droplets(
      method = "supplied",
      is_empty_column = "is_empty"
    ),
    .verbose = FALSE
  ),
  pattern = "must be a logical column",
  info = "resolve empty droplets: non-logical column errors"
)

obs_with_na <- data.table::copy(obs_resolve)
obs_with_na[1L, is_empty := NA]

expect_error(
  current = bixverse:::.resolve_empty_droplets(
    obs = obs_with_na,
    empty_params = params_sc_empty_droplets(
      method = "supplied",
      is_empty_column = "is_empty"
    ),
    .verbose = FALSE
  ),
  pattern = "contains NAs",
  info = "resolve empty droplets: NAs in the mask error"
)

## parameter constructors ------------------------------------------------------

expect_true(
  current = checkmate::testClass(
    params_sc_cellsweep(),
    "params_sc_cellsweep"
  ),
  info = "cellsweep params: constructor tags the class"
)

expect_true(
  current = isTRUE(bixverse:::checkScCellsweep(params_sc_cellsweep())),
  info = "cellsweep params: the default list passes its own check"
)

expect_true(
  current = is.character(
    bixverse:::checkScCellsweep(list(freeze_empties = TRUE))
  ),
  info = "cellsweep params: an incomplete list is rejected"
)

expect_error(
  current = params_sc_cellsweep(freeze_empties = FALSE),
  pattern = "is not supported",
  info = "cellsweep params: unfreezing the empty droplets is refused"
)

expect_true(
  current = isTRUE(
    bixverse:::checkScEmptyDroplets(
      params_sc_empty_droplets(
        method = "supplied",
        is_empty_column = "is_empty"
      )
    )
  ),
  info = "empty droplet params: a valid list passes its own check"
)

expect_true(
  current = is.character(
    bixverse:::checkScEmptyDroplets(list(method = "knee"))
  ),
  info = "empty droplet params: an incomplete list is rejected"
)

expect_error(
  current = params_sc_empty_droplets(method = "supplied"),
  pattern = "needs `is_empty_column`",
  info = "empty droplet params: supplied without its column errors"
)

expect_error(
  current = params_sc_empty_droplets(method = "umi_cutoff"),
  pattern = "needs `umi_cutoff`",
  info = "empty droplet params: umi cutoff without its cutoff errors"
)

expect_error(
  current = params_sc_empty_droplets(method = "expected_cells"),
  pattern = "needs `expected_cells`",
  info = "empty droplet params: expected cells without its count errors"
)

## synthetic data generator ----------------------------------------------------

expect_true(
  current = checkmate::testClass(fixture$counts, "dgRMatrix"),
  info = "cellsweep fixture: counts come back as a dgRMatrix"
)

expect_true(
  current = sum(is_empty_true) == 2000L && nrow(fixture$counts) == 2600L,
  info = "cellsweep fixture: real barcodes first, empty droplets after"
)

expect_true(
  current = all(is.na(fixture$obs$cell_grp[is_empty_true])),
  info = "cellsweep fixture: empty droplets carry no cell type label"
)

expect_true(
  current = abs(sum(fixture$ambient_true) - 1) < 1e-9,
  info = "cellsweep fixture: the planted soup is a distribution"
)

expect_true(
  current = all(
    abs(rowSums(fixture$celltype_profiles_true) - 1) < 1e-9
  ),
  info = "cellsweep fixture: every cell type profile is a distribution"
)

expect_error(
  current = generate_cellsweep_test_data(
    syn_data_params = params_sc_synthetic_cellsweep(n_empty = 40L, n_real = 0L)
  ),
  info = "cellsweep fixture: an empty experiment is refused"
)

## the s7 method ---------------------------------------------------------------

sc_input <- cellsweep_test_object(test_dir_input, fixture)

sc_denoised <- cellsweep_sc(
  target = SingleCells(dir_data = test_dir_output),
  input = sc_input,
  celltype_column = "cell_grp",
  sample_column = "sample_id",
  empty_params = params_sc_empty_droplets(
    method = "supplied",
    is_empty_column = "is_empty"
  ),
  streaming = 0L,
  .verbose = FALSE
)

obs_denoised <- get_sc_obs(sc_denoised, filtered = FALSE)
var_denoised <- get_sc_var(sc_denoised)

expect_true(
  current = all(S7::prop(sc_denoised, "dims") == c(n_real, n_genes)),
  info = "cellsweep: only the annotated barcodes make it into the output"
)

expect_true(
  current = checkmate::testNames(
    names(obs_denoised),
    must.include = c(
      "cellsweep_alpha",
      "cellsweep_z",
      "cellsweep_beta",
      "cellsweep_ll",
      "cellsweep_converged"
    )
  ),
  info = "cellsweep: the fit lands in the obs table"
)

expect_true(
  current = checkmate::testNames(
    names(var_denoised),
    must.include = c("cellsweep_ambient", "no_cells_exp")
  ),
  info = "cellsweep: the ambient profile lands in the var table"
)

expect_true(
  current = all(obs_denoised$cellsweep_converged),
  info = "cellsweep: the EM converged"
)

expect_true(
  current = checkmate::qtest(
    obs_denoised$cellsweep_alpha,
    sprintf("N%s[0,1]", n_real)
  ),
  info = "cellsweep: one ambient fraction per barcode, inside [0, 1]"
)

expect_true(
  current = cor(
    obs_denoised$cellsweep_alpha,
    obs_denoised$alpha_true,
    method = "spearman"
  ) >=
    0.9,
  info = "cellsweep: the fitted ambient fraction tracks the planted one"
)

expect_true(
  current = abs(
    mean(obs_denoised$cellsweep_alpha - obs_denoised$alpha_true)
  ) <
    0.05,
  info = "cellsweep: the ambient fraction is not systematically off"
)

expect_true(
  current = cor(var_denoised$cellsweep_ambient, fixture$ambient_true) >= 0.99,
  info = "cellsweep: the estimated soup matches the planted one"
)

expect_true(
  current = mean(obs_denoised$cellsweep_z == obs_denoised$cell_grp) == 1,
  info = "cellsweep: no barcode gets reassigned to another cell type"
)

### denoising ------------------------------------------------------------------

raw_type_2 <- fixture$counts[which(fixture$obs$cell_grp == "cell_type_2"), ]
clean_type_2 <- get_sc_counts(
  sc_denoised,
  assay = "raw",
  cell_indices = which(obs_denoised$cell_grp == "cell_type_2"),
  use_cells_to_keep = FALSE,
  .verbose = FALSE
)

expect_true(
  current = soup_fraction(clean_type_2) < 0.25 * soup_fraction(raw_type_2),
  info = "cellsweep: the soup is largely gone from the wrong cell type"
)

expect_true(
  current = sum(clean_type_2) < sum(raw_type_2),
  info = "cellsweep: denoising removes counts rather than adding them"
)

expect_true(
  current = all(clean_type_2 >= 0),
  info = "cellsweep: the denoised counts stay non-negative"
)

### inferring the mask inside the method ---------------------------------------

sc_inferred <- cellsweep_sc(
  target = SingleCells(dir_data = test_dir_inferred),
  input = sc_input,
  celltype_column = "cell_grp",
  sample_column = "sample_id",
  empty_params = params_sc_empty_droplets(
    method = "umi_cutoff",
    umi_cutoff = 500L
  ),
  streaming = 0L,
  .verbose = FALSE
)

expect_equivalent(
  current = get_sc_obs(sc_inferred, filtered = FALSE)$cellsweep_alpha,
  target = obs_denoised$cellsweep_alpha,
  info = "cellsweep: inferring the mask reproduces the supplied one"
)

### error handling -------------------------------------------------------------

expect_error(
  current = cellsweep_sc(
    target = sc_input,
    input = sc_input,
    celltype_column = "cell_grp",
    sample_column = "sample_id",
    empty_params = params_sc_empty_droplets(
      method = "supplied",
      is_empty_column = "is_empty"
    ),
    .verbose = FALSE
  ),
  pattern = "must point at different directories",
  info = "cellsweep: writing back into the input directory errors"
)

expect_error(
  current = cellsweep_sc(
    target = SingleCells(dir_data = test_dir_output),
    input = sc_input,
    celltype_column = "not_a_column",
    sample_column = "sample_id",
    empty_params = params_sc_empty_droplets(
      method = "supplied",
      is_empty_column = "is_empty"
    ),
    .verbose = FALSE
  ),
  pattern = "is not a column in the obs table",
  info = "cellsweep: a missing cell type column errors"
)

expect_error(
  current = cellsweep_sc(
    target = SingleCells(dir_data = test_dir_output),
    input = sc_input,
    celltype_column = "cell_grp",
    sample_column = "not_a_column",
    empty_params = params_sc_empty_droplets(
      method = "supplied",
      is_empty_column = "is_empty"
    ),
    .verbose = FALSE
  ),
  pattern = "is not a column in the obs table",
  info = "cellsweep: a missing sample column errors"
)

obs_unlabelled <- data.table::copy(fixture$obs)
obs_unlabelled[, cell_grp := NA_character_]

sc_unlabelled <- cellsweep_test_object(
  sc_test_dir("cellsweep_unlabelled"),
  fixture,
  obs = obs_unlabelled
)

expect_error(
  current = cellsweep_sc(
    target = SingleCells(dir_data = test_dir_output),
    input = sc_unlabelled,
    celltype_column = "cell_grp",
    sample_column = "sample_id",
    empty_params = params_sc_empty_droplets(
      method = "supplied",
      is_empty_column = "is_empty"
    ),
    .verbose = FALSE
  ),
  pattern = "No barcode is both annotated and passing QC",
  info = "cellsweep: an object without labels errors"
)

# the load-time cutoffs are irreversible, so an object ingested with real QC
# has no empty droplets left to fit the ambient profile on
sc_filtered <- load_r_data(
  object = SingleCells(dir_data = test_dir_filtered),
  counts = fixture$counts,
  obs = fixture$obs,
  var = fixture$var,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 0L,
    min_lib_size = 500L,
    min_cells = 0L
  ),
  streaming = 0L,
  .verbose = FALSE
)

expect_error(
  current = cellsweep_sc(
    target = SingleCells(dir_data = test_dir_output),
    input = sc_filtered,
    celltype_column = "cell_grp",
    sample_column = "sample_id",
    empty_params = params_sc_empty_droplets(
      method = "supplied",
      is_empty_column = "is_empty"
    ),
    .verbose = FALSE
  ),
  pattern = "the empty droplets were filtered out at ingest",
  info = "cellsweep: an object ingested with QC errors"
)

# clean up ---------------------------------------------------------------------

sc_test_cleanup(
  test_dir_input,
  test_dir_output,
  test_dir_inferred,
  test_dir_filtered,
  file.path(tempdir(), "cellsweep_unlabelled")
)
