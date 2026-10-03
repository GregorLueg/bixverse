# h5ad export ------------------------------------------------------------------

source("helper_sc.R", local = TRUE)

test_temp_dir <- sc_test_dir("io_h5_write")
roundtrip_dir <- sc_test_dir("io_h5_write_roundtrip")

## fixture ---------------------------------------------------------------------

fixture <- sc_test_fixture()

# a repeating character column exercises the categorical encoding, a unique one
# the string array path
obs <- data.table::copy(fixture$obs)
obs[, "unique_tag" := sprintf("tag_%04d", seq_len(.N))]
obs[, "is_even" := seq_len(.N) %% 2L == 0L]
# missing values in a categorical and in an integer column, both legal h5ad
obs[, "grp_na" := c("a", "b", NA)[seq_len(.N) %% 3L + 1L]]
obs[, "n_na" := ifelse(seq_len(.N) %% 4L == 0L, NA_integer_, seq_len(.N))]

sc_object <- sc_test_object(dir = test_temp_dir, fixture = fixture, obs = obs)
sc_object <- sc_test_prepped(object = sc_object, fixture = fixture)

counts_expected <- sc_object[,, return_format = "cell"]
obs_expected <- sc_object[[]]
var_expected <- get_sc_var(sc_object)

h5_out <- file.path(test_temp_dir, "export.h5ad")

res_path <- save_h5ad(
  object = sc_object,
  h5_path = h5_out,
  chunk_size = 137L,
  .verbose = FALSE
)

# tests ------------------------------------------------------------------------

## file creation ---------------------------------------------------------------

expect_true(
  current = file.exists(h5_out),
  info = "h5ad export - file is written"
)

expect_equal(
  current = res_path,
  target = path.expand(h5_out),
  info = "h5ad export - the file path is returned"
)

expect_error(
  current = save_h5ad(
    object = sc_object,
    h5_path = h5_out,
    overwrite = FALSE,
    .verbose = FALSE
  ),
  info = "h5ad export - refuses to overwrite when told not to"
)

## anndata encoding ------------------------------------------------------------

# anndata refuses to read a file whose groups are not tagged, so the encoding
# attributes are the actual contract with ScanPy

root_attrs <- rhdf5::h5readAttributes(h5_out, "/")

expect_equal(
  current = root_attrs[["encoding-type"]],
  target = "anndata",
  info = "h5ad export - the root carries the anndata encoding"
)

x_attrs <- rhdf5::h5readAttributes(h5_out, "X")

expect_equal(
  current = x_attrs[["encoding-type"]],
  target = "csr_matrix",
  info = "h5ad export - X is tagged as a CSR matrix"
)

expect_equal(
  current = as.integer(x_attrs[["shape"]]),
  target = dim(counts_expected),
  info = "h5ad export - the X shape attribute matches the object"
)

obs_attrs <- rhdf5::h5readAttributes(h5_out, "obs")

expect_equal(
  current = obs_attrs[["encoding-type"]],
  target = "dataframe",
  info = "h5ad export - obs is tagged as a dataframe"
)

expect_equal(
  current = as.vector(rhdf5::h5read(h5_out, "obs/_index")),
  target = obs_expected$cell_id,
  info = "h5ad export - the cell ids are the obs index"
)

expect_true(
  current = setequal(
    obs_attrs[["column-order"]],
    setdiff(names(obs_expected), c("cell_id", "cell_idx", "to_keep"))
  ),
  info = "h5ad export - column-order lists every obs column but the index"
)

expect_false(
  current = "to_keep" %in% obs_attrs[["column-order"]],
  info = "h5ad export - bixverse bookkeeping columns are not exported"
)

expect_equal(
  current = as.vector(rhdf5::h5read(h5_out, "var/_index")),
  target = var_expected$gene_id,
  info = "h5ad export - the gene ids are the var index"
)

for (grp in c("obsm", "varm", "obsp", "uns")) {
  expect_equal(
    current = rhdf5::h5readAttributes(h5_out, grp)[["encoding-type"]],
    target = "dict",
    info = sprintf("h5ad export - %s is tagged as a dict", grp)
  )
}

## column encodings ------------------------------------------------------------

# a repeating character column becomes a categorical, a column of unique
# strings stays a string array

expect_equal(
  current = rhdf5::h5readAttributes(h5_out, "obs/cell_grp")[["encoding-type"]],
  target = "categorical",
  info = "h5ad export - a repeating character column becomes categorical"
)

expect_equal(
  current = rhdf5::h5readAttributes(
    h5_out,
    "obs/unique_tag"
  )[["encoding-type"]],
  target = "string-array",
  info = "h5ad export - a column of unique strings stays a string array"
)

expect_equal(
  current = rhdf5::h5readAttributes(
    h5_out,
    "obs/cell_grp/codes"
  )[["encoding-type"]],
  target = "array",
  info = "h5ad export - the categorical codes are tagged as an array"
)

cat_levels <- as.vector(rhdf5::h5read(h5_out, "obs/cell_grp/categories"))
cat_codes <- as.vector(rhdf5::h5read(h5_out, "obs/cell_grp/codes"))

expect_equal(
  current = cat_levels[cat_codes + 1L],
  target = obs_expected$cell_grp,
  info = "h5ad export - categorical codes and categories rebuild the column"
)

expect_equal(
  current = as.vector(rhdf5::h5read(h5_out, "obs/unique_tag")),
  target = obs_expected$unique_tag,
  info = "h5ad export - string array values survive"
)

# rhdf5 hands an HDF5 boolean enum back as a factor over "FALSE" / "TRUE";
# h5py and therefore ScanPy see a numpy bool
expect_equal(
  current = as.logical(as.character(rhdf5::h5read(h5_out, "obs/is_even"))),
  target = obs_expected$is_even,
  info = "h5ad export - logical columns survive the boolean enum"
)

na_codes <- as.vector(rhdf5::h5read(h5_out, "obs/grp_na/codes"))

expect_equal(
  current = na_codes == -1L,
  target = is.na(obs_expected$grp_na),
  info = "h5ad export - a missing categorical value is written as code -1"
)

n_na_written <- as.vector(rhdf5::h5read(h5_out, "obs/n_na"))

# NA_integer_ has no float representation but NaN; rhdf5 reads NaN back as NA
expect_true(
  current = is.double(n_na_written),
  info = "h5ad export - an integer column with NA is widened to float"
)

expect_equal(
  current = is.na(n_na_written),
  target = is.na(obs_expected$n_na),
  info = "h5ad export - a missing integer stays missing"
)

## counts ----------------------------------------------------------------------

indptr <- as.vector(rhdf5::h5read(h5_out, "X/indptr"))

expect_equal(
  current = length(indptr),
  target = nrow(counts_expected) + 1L,
  info = "h5ad export - one index pointer per cell plus the terminator"
)

expect_equal(
  current = indptr,
  target = as.numeric(counts_expected@p),
  info = "h5ad export - the streamed index pointer matches the object"
)

expect_equal(
  current = as.vector(rhdf5::h5read(h5_out, "X/indices")),
  target = counts_expected@j,
  info = "h5ad export - the streamed indices match the object"
)

expect_equal(
  current = as.vector(rhdf5::h5read(h5_out, "X/data")),
  target = counts_expected@x,
  info = "h5ad export - the streamed data match the object"
)

### chunking invariance --------------------------------------------------------

# the streaming pointer arithmetic is the easy thing to get wrong, so the
# result has to be independent of how the cells are batched

h5_single_chunk <- file.path(test_temp_dir, "export_single_chunk.h5ad")

save_h5ad(
  object = sc_object,
  h5_path = h5_single_chunk,
  chunk_size = 1000000L,
  .verbose = FALSE
)

for (slot in c("X/data", "X/indices", "X/indptr")) {
  expect_equal(
    current = as.vector(rhdf5::h5read(h5_out, slot)),
    target = as.vector(rhdf5::h5read(h5_single_chunk, slot)),
    info = sprintf("h5ad export - %s does not depend on the chunk size", slot)
  )
}

## cached artefacts ------------------------------------------------------------

pca_expected <- get_embedding(sc_object, embd_name = "pca")

# rhdf5 hands back the transpose of the on-disk (obs x k) layout
pca_written <- t(rhdf5::h5read(h5_out, "obsm/X_pca"))

expect_equal(
  current = pca_written,
  target = unname(pca_expected),
  info = "h5ad export - the PCA factors land in obsm/X_pca"
)

expect_equal(
  current = dim(t(rhdf5::h5read(h5_out, "varm/PCs"))),
  target = c(nrow(var_expected), ncol(pca_expected)),
  info = "h5ad export - the PCA loadings are padded out to every gene"
)

loadings_written <- t(rhdf5::h5read(h5_out, "varm/PCs"))
hvg_idx <- get_hvg(sc_object) + 1L

expect_equal(
  current = loadings_written[hvg_idx, ],
  target = unname(get_pca_loadings(sc_object)),
  info = "h5ad export - the loadings sit on the HVG rows"
)

expect_true(
  current = all(loadings_written[-hvg_idx, ] == 0),
  info = "h5ad export - non-HVG loading rows are zero"
)

expect_equal(
  current = rhdf5::h5readAttributes(
    h5_out,
    "obsp/connectivities"
  )[["encoding-type"]],
  target = "csr_matrix",
  info = "h5ad export - the sNN graph lands in obsp as a CSR matrix"
)

expect_equal(
  current = as.integer(
    rhdf5::h5readAttributes(h5_out, "obsp/connectivities")[["shape"]]
  ),
  target = rep(nrow(counts_expected), 2L),
  info = "h5ad export - the sNN graph is square over the cells"
)

expect_equal(
  current = as.vector(rhdf5::h5read(
    h5_out,
    "uns/neighbors/connectivities_key"
  )),
  target = "connectivities",
  info = "h5ad export - ScanPy is pointed at the connectivities"
)

# anndata reads a scalar string as `elem.asstr()[()]`, so a rank one dataset
# would come back as a length one array rather than as the string itself

key_fid <- rhdf5::H5Fopen(h5_out)
key_did <- rhdf5::H5Dopen(key_fid, "uns/neighbors/connectivities_key")
key_sid <- rhdf5::H5Dget_space(key_did)
key_rank <- rhdf5::H5Sget_simple_extent_dims(key_sid)$rank
rhdf5::H5Sclose(key_sid)
rhdf5::H5Dclose(key_did)
rhdf5::H5Fclose(key_fid)

expect_equal(
  current = key_rank,
  target = 0L,
  info = "h5ad export - a scalar string is a rank zero dataset"
)

## round trip ------------------------------------------------------------------

# the strongest available check: read the exported file back through the
# package's own h5ad ingestion and compare

sc_roundtrip <- SingleCells(dir_data = roundtrip_dir)

sc_roundtrip <- stream_h5ad(
  object = sc_roundtrip,
  h5_path = h5_out,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 0L,
    min_lib_size = 0L,
    min_cells = 0L
  ),
  .verbose = FALSE
)

counts_roundtrip <- sc_roundtrip[,, return_format = "cell"]

expect_equal(
  current = dim(counts_roundtrip),
  target = dim(counts_expected),
  info = "h5ad round trip - the dimensions survive"
)

expect_equal(
  current = as.vector(counts_roundtrip),
  target = as.vector(counts_expected),
  info = "h5ad round trip - the counts survive"
)

expect_equal(
  current = get_cell_names(sc_roundtrip),
  target = obs_expected$cell_id,
  info = "h5ad round trip - the cell identifiers survive"
)

expect_equal(
  current = get_gene_names(sc_roundtrip),
  target = var_expected$gene_id,
  info = "h5ad round trip - the gene identifiers survive"
)

expect_equal(
  current = as.character(sc_roundtrip[[]]$cell_grp),
  target = as.character(obs_expected$cell_grp),
  info = "h5ad round trip - a categorical obs column survives"
)

## cell filtering --------------------------------------------------------------

# the export follows `cells_to_keep`, i.e. it matches what `object[]` returns

keep <- seq_len(200L)

sc_filtered <- suppressWarnings(
  set_cells_to_keep(sc_object, cells_to_keep = as.integer(keep))
)

h5_filtered <- file.path(test_temp_dir, "export_filtered.h5ad")

suppressWarnings(
  save_h5ad(
    object = sc_filtered,
    h5_path = h5_filtered,
    .verbose = FALSE
  )
)

expect_equal(
  current = as.integer(
    rhdf5::h5readAttributes(h5_filtered, "X")[["shape"]]
  ),
  target = c(length(keep), ncol(counts_expected)),
  info = "h5ad export - only the kept cells are written"
)

expect_equal(
  current = as.vector(rhdf5::h5read(h5_filtered, "obs/_index")),
  target = obs_expected$cell_id[keep],
  info = "h5ad export - the obs table follows the kept cells"
)

counts_kept <- methods::as(
  methods::as(counts_expected[keep, ], "RsparseMatrix"),
  "dgRMatrix"
)

expect_equal(
  current = as.vector(rhdf5::h5read(h5_filtered, "X/data")),
  target = counts_kept@x,
  info = "h5ad export - the counts follow the kept cells"
)

# the cached PCA covers every cell, so it is stale for the filtered object and
# must be left out rather than written at the wrong length
expect_true(
  current = !("X_pca" %in%
    rhdf5::h5ls(h5_filtered)[
      rhdf5::h5ls(h5_filtered)$group == "/obsm",
      "name"
    ]),
  info = "h5ad export - a stale embedding is skipped, not written"
)

## column types --------------------------------------------------------------

# which anndata encoding each obs/var column gets is decided in R; the Rust
# side only maps the types across

cols <- bixverse:::.h5ad_columns(data.table::data.table(
  grp = c("a", "b", "a"),
  tag = c("x", "y", "z"),
  all_na = NA_character_,
  int_na = c(1L, NA, 3L),
  lgl_na = c(TRUE, NA, FALSE),
  lgl = c(TRUE, FALSE, TRUE),
  day = as.Date("2026-01-01") + 0:2
))

expect_true(
  current = is.factor(cols$grp),
  info = "h5ad columns - a repeating character column becomes a factor"
)

expect_true(
  current = is.character(cols$tag),
  info = "h5ad columns - unique strings stay a string array"
)

expect_equal(
  current = cols$all_na,
  target = rep(NA_character_, 3L),
  info = "h5ad columns - an all NA column has no category and stays strings"
)

expect_true(
  current = is.double(cols$int_na) && is.na(cols$int_na[2L]),
  info = "h5ad columns - integers with NA are widened to double"
)

expect_true(
  current = is.double(cols$lgl_na),
  info = "h5ad columns - logicals with NA are widened to double"
)

expect_true(
  current = is.logical(cols$lgl),
  info = "h5ad columns - complete logicals stay booleans"
)

expect_equal(
  current = cols$day,
  target = c("2026-01-01", "2026-01-02", "2026-01-03"),
  info = "h5ad columns - unsupported types are stringified, not dropped"
)

## input checks --------------------------------------------------------------

# a bad input has to fail before the file exists, not after the counts were
# streamed into it

h5_bad <- file.path(test_temp_dir, "bad.h5ad")

expect_error(
  current = rs_save_h5ad(
    f_path_cells = bixverse:::get_rust_count_cell_f_path(sc_object),
    h5_path = h5_bad,
    cell_indices = 0:9,
    norm = FALSE,
    obs_index = sprintf("c%i", 1:10),
    obs = list(too_short = 1:5),
    var_index = var_expected$gene_id,
    var = list(),
    obsm = list(),
    varm = list(),
    obsp = list(),
    uns_json = "{}",
    chunk_size = 5L
  ),
  pattern = "too_short",
  info = "h5ad export - a column of the wrong length is rejected"
)

expect_false(
  current = file.exists(h5_bad),
  info = "h5ad export - a rejected input leaves no file behind"
)

# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir, roundtrip_dir)
