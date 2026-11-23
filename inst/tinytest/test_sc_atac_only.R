# atac only experiments ---------------------------------------------------------

# port of bixverse-rs python/tests/test_atac_only.py: SingleCellsAtac shares its
# implementation with the multi-modal class, so most checks are that both give
# the same numbers on the same cells

source("helper_atac.R", local = TRUE)
if (!ATAC_OK) {
  exit_file("ATAC fixtures not found; set BIXVERSE_ATAC_FIXTURES")
}

## load ------------------------------------------------------------------------

x <- atac_only_experiment()
mm <- atac_mm_experiment()
signac <- atac_tsv("signac/cell_metrics.tsv")

obs <- get_sc_obs(x, cols = c("cell_id", "lib_size"))
expect_equal(obs$cell_id, signac$barcode, info = "barcode order")
expect_equal(
  as.integer(obs$lib_size),
  signac$nCount_peaks,
  info = "counts match Signac"
)
expect_equal(x@atac_info$kind, "peak", info = "peak store")
expect_equal(x@atac_info$count_mode, "fragment", info = "count mode")
expect_false(x@atac_info$tf_idf, info = "no TF-IDF yet")
ours <- get_sc_var(x, modality = "atac")
expect_equal(
  ours$gene_id,
  get_sc_var(mm, modality = "atac")$feature_id,
  info = "same features as the multi-modal store"
)

## barcode sources -------------------------------------------------------------

counts <- barcode_fragment_counts(ATAC_FRAGS)
expect_false(
  is.unsorted(rev(counts$n_fragments)),
  info = "most fragments first"
)
cut <- counts$n_fragments[101]
dir_cut <- tempfile("atac_cut")
dir.create(dir_cut)
x_cut <- load_atac_fragments(
  SingleCellsAtac(dir_cut),
  ATAC_FRAGS,
  peaks = ATAC_PEAKS,
  min_fragments = cut,
  .verbose = FALSE
)
expect_equal(
  get_cell_names(x_cut),
  sort(counts$barcode[counts$n_fragments >= cut], method = "radix"),
  info = "min_fragments calls the cells"
)

csv <- tempfile(fileext = ".csv")
writeLines(
  c(
    "barcode,is__cell_barcode",
    "NO_BARCODE,0",
    sprintf(
      "%s,%i",
      ATAC_BARCODES,
      as.integer((seq_along(ATAC_BARCODES) - 1L) %% 2L == 0L)
    )
  ),
  csv
)
expect_equal(
  read_cellranger_barcodes(csv),
  ATAC_BARCODES[seq(1L, length(ATAC_BARCODES), by = 2L)],
  info = "Cell Ranger singlecell.csv"
)

dir_bad <- tempfile("atac_bad")
dir.create(dir_bad)
expect_error(
  load_atac_fragments(SingleCellsAtac(dir_bad), ATAC_FRAGS, peaks = ATAC_PEAKS),
  pattern = "exactly one",
  info = "no cell source"
)
expect_error(
  load_atac_fragments(
    SingleCellsAtac(dir_bad),
    ATAC_FRAGS,
    peaks = ATAC_PEAKS,
    barcodes = ATAC_BARCODES,
    min_fragments = 10L
  ),
  pattern = "exactly one",
  info = "two cell sources"
)

## barcode counts against fragment qc ------------------------------------------

xq <- fragment_qc_sc(x, ATAC_FRAGS, ATAC_GTF, mito_chroms = character())
counts_all <- barcode_fragment_counts(ATAC_FRAGS, exclude_chroms = character())
want <- counts_all$n_fragments[match(get_cell_names(x), counts_all$barcode)]
want[is.na(want)] <- 0
expect_equal(
  get_sc_obs(xq)$atac_n_fragments,
  as.numeric(want),
  info = "barcode counts agree with the QC pass"
)

## one lsi code path -----------------------------------------------------------

lsi_of <- function(e) {
  e <- run_tfidf_sc(e)
  e <- find_top_features_sc(e, min_cutoff = "q5")
  calculate_lsi_sc(
    e,
    lsi_params = params_atac_lsi(n_comp = 10L, solver = "exact"),
    .verbose = FALSE
  )
}
x_lsi <- lsi_of(x)
mm_lsi <- lsi_of(mm)
expect_equal(
  unname(get_embedding(x_lsi, "lsi")),
  unname(get_embedding(mm_lsi, "lsi", modality = "atac")),
  tolerance = 1e-5,
  info = "same LSI on both classes"
)
expect_equal(
  x_lsi@atac_info$depth_cor,
  mm_lsi@atac_info$depth_cor,
  tolerance = 1e-6,
  info = "same depth correlation"
)

## modality names the root -----------------------------------------------------

x_nn <- find_neighbours_sc(
  x_lsi,
  embd_to_use = "lsi",
  modality = "atac",
  .verbose = FALSE
)
expect_equal(nrow(get_knn_mat(x_nn)), 500L, info = "atac is the root cache")
expect_error(
  find_neighbours_sc(
    x_lsi,
    embd_to_use = "lsi",
    modality = "adt",
    .verbose = FALSE
  ),
  pattern = "no 'adt' modality",
  info = "no ADT on an ATAC-only object"
)

## differential accessibility --------------------------------------------------

a <- ATAC_BARCODES[1:250]
b <- ATAC_BARCODES[251:500]
ours <- find_markers_sc(x, a, b, .verbose = FALSE)
theirs <- find_markers_sc(mm, a, b, modality = "atac", .verbose = FALSE)
expect_equal(names(ours)[1], "feature_id", info = "keyed by feature id")
expect_equal(ours, theirs, info = "LR by default, same on both classes")

## rna methods refuse peaks ----------------------------------------------------

expect_error(find_hvg_sc(x), pattern = "works on genes", info = "no HVG")
expect_error(
  calculate_pca_sc(x, no_pcs = 5L),
  pattern = "works on genes",
  info = "no PCA"
)
expect_error(scrublet_sc(x), pattern = "needs RNA counts", info = "no scrublet")

## peak accessibility plots ----------------------------------------------------

peak <- get_sc_var(x, modality = "atac")$gene_id[1]
for (e in list(x, mm)) {
  expr <- extract_gene_expression(e, features = peak, modality = "atac")
  expect_equal(names(expr), c("cell_id", peak), info = "one column per peak")
  expect_equal(nrow(expr), 500L, info = "every kept cell")
}

## state is rebuilt from the store ---------------------------------------------

save_sc_exp_to_disk(x_lsi)
back <- load_existing(SingleCellsAtac(x_lsi@dir_data), .verbose = FALSE)
expect_true(back@atac_info$tf_idf, info = "TF-IDF from the header")
expect_equal(back@atac_info$count_mode, "fragment", info = "count mode")
expect_equal(
  back@atac_info$depth_cor,
  x_lsi@atac_info$depth_cor,
  info = "LSI fit from the memory file"
)
expect_true(
  "lsi" %in% get_available_embeddings(back),
  info = "LSI embedding restored"
)

unlink(file.path(x_lsi@dir_data, "memory.qs2"))
bare_back <- load_existing(SingleCellsAtac(x_lsi@dir_data), .verbose = FALSE)
expect_true(bare_back@atac_info$tf_idf, info = "TF-IDF without a memory file")
mm_back <- load_existing(
  SingleCellsMultiModal(mm@dir_data),
  .verbose = FALSE
)
expect_equal(
  mm_back@atac_info$count_mode,
  "fragment",
  info = "multi-modal ATAC state from the store"
)
expect_error(
  load_existing(SingleCellsMultiModal(x@dir_data), .verbose = FALSE),
  pattern = "SingleCellsAtac",
  info = "an ATAC-only store is not multi-modal"
)

## gene activity follows the filter --------------------------------------------

xg <- atac_only_experiment()
ga <- gene_activity_sc(xg, ATAC_FRAGS, ATAC_GTF, .verbose = FALSE)
expect_equal(length(get_cells_to_keep(ga)), 500L, info = "all cells")
xg <- set_cells_to_keep(xg, ATAC_BARCODES[1:100])
expect_equal(
  length(get_cells_to_keep(get_gene_activity(xg))),
  100L,
  info = "get_gene_activity() takes the filter"
)
ga_db <- bixverse:::SingleCellDuckDB$new(
  db_dir = file.path(xg@dir_data, "gene_activity")
)
expect_equal(
  length(ga_db$get_cells_to_keep()),
  100L,
  info = "the filter reaches the gene activity store"
)

## frip after a recount --------------------------------------------------------

mf <- fragment_qc_sc(mm, ATAC_FRAGS, ATAC_GTF)
before <- get_sc_obs(mf)$atac_frip
mf <- add_atac_fragments_sc(
  mf,
  ATAC_FRAGS,
  peaks = ATAC_PEAKS,
  mode = "insertion",
  overwrite = TRUE,
  .verbose = FALSE
)
obs <- get_sc_obs(mf)
ok <- obs$atac_n_fragments > 0
expect_equal(
  obs$atac_frip[ok],
  obs$atac_lib_size[ok] / (2 * obs$atac_n_fragments[ok]),
  info = "FRiP recomputed against the new counts"
)
expect_false(
  isTRUE(all.equal(before[ok], obs$atac_frip[ok])),
  info = "FRiP moved"
)

## mismatched inputs fail loudly -----------------------------------------------

dir_sfx <- tempfile("atac_sfx")
dir.create(dir_sfx)
expect_error(
  load_atac_fragments(
    SingleCellsAtac(dir_sfx),
    ATAC_FRAGS,
    peaks = ATAC_PEAKS,
    barcodes = sub("-1$", "", ATAC_BARCODES),
    .verbose = FALSE
  ),
  pattern = "barcodes are not in the input",
  info = "no barcode in the file"
)
plain <- get_sc_obs(fragment_qc_sc(x, ATAC_FRAGS, ATAC_GTF))
mapped <- get_sc_obs(fragment_qc_sc(
  x,
  ATAC_FRAGS,
  ATAC_GTF,
  barcode_map = stats::setNames(ATAC_BARCODES, ATAC_BARCODES)
))
expect_equal(plain, mapped, info = "identity barcode map")
expect_error(
  fragment_qc_sc(
    x,
    ATAC_FRAGS,
    ATAC_GTF,
    barcode_map = \(ids) sub("-1$", "", ids)
  ),
  pattern = "barcodes are not in the input",
  info = "a wrong barcode map"
)
dir_chr <- tempfile("atac_chr")
dir.create(dir_chr)
expect_error(
  load_atac_fragments(
    SingleCellsAtac(dir_chr),
    ATAC_FRAGS,
    peaks = data.table::data.table(
      chrom = "22",
      start = 16000000L,
      end = 16001000L
    ),
    barcodes = ATAC_BARCODES,
    .verbose = FALSE
  ),
  pattern = "chromosome naming",
  info = "chromosome names do not match"
)

## profiles, annotation and blacklist ------------------------------------------

prof <- fragment_profiles_sc(x, ATAC_FRAGS, ATAC_GTF)
expect_equal(nrow(prof$tss), 4001L, info = "TSS profile width")
expect_equal(nrow(prof$lengths), 1000L, info = "fragment length bins")
expect_true(
  mean(prof$tss$enrichment[abs(prof$tss$offset) < 100]) > 3,
  info = "TSS enrichment at the centre"
)
expect_equal(
  sum(prof$lengths$fraction),
  1,
  tolerance = 1e-9,
  info = "fractions"
)

xa <- annotate_peaks_sc(x, ATAC_GTF)
var <- get_sc_var(xa)
expect_true("closest_gene" %in% names(var), info = "closest gene in var")
expect_true(all(var$distance_to_gene >= 0), info = "distances")

peaks_dt <- data.table::fread(
  ATAC_PEAKS,
  header = FALSE,
  col.names = c("chrom", "start", "end")
)
kept <- remove_blacklisted(peaks_dt, file.path(ATAC_FIX, "blacklist.bed"))
expect_true(
  nrow(kept) > 0 && nrow(kept) <= nrow(peaks_dt),
  info = "blacklist filter"
)

## merge fragments -------------------------------------------------------------

out <- tempfile(fileext = ".tsv.gz")
summary <- merge_fragments(
  c(ATAC_FRAGS, ATAC_FRAGS),
  out,
  prefixes = c("s1_", "s2_")
)
expect_true(
  summary$n_fragments[1] == summary$n_fragments[2] &&
    summary$n_fragments[1] > 0,
  info = "both inputs merged"
)
expect_true(file.exists(paste0(out, ".tbi")), info = "tabix index written")
dir_m <- tempfile("atac_merged")
dir.create(dir_m)
merged <- load_atac_fragments(
  SingleCellsAtac(dir_m),
  out,
  peaks = ATAC_PEAKS,
  barcodes = c(
    paste0("s1_", ATAC_BARCODES[1:50]),
    paste0("s2_", ATAC_BARCODES[1:50])
  ),
  .verbose = FALSE
)
lib <- get_sc_obs(merged)$lib_size
expect_equal(lib[1:50], lib[51:100], info = "samples stay apart")

## peak to gene ----------------------------------------------------------------

get_sc_duckdb(mm)$add_data_obs(data.table::data.table(
  grp = sprintf("g%i", (seq_len(500L) - 1L) %% 20L)
))
links <- peak_to_gene_sc(mm, ATAC_GTF, "grp", fdr = NULL, .verbose = FALSE)
expect_equal(
  sort(names(links)),
  sort(c("peak_id", "gene_id", "gene_name", "distance", "r", "p", "fdr")),
  info = "peak to gene columns"
)
filtered <- peak_to_gene_sc(mm, ATAC_GTF, "grp", fdr = 0.5, .verbose = FALSE)
expect_equal(
  nrow(filtered),
  sum(links$fdr < 0.5),
  info = "the FDR filter keeps the links below the cutoff"
)
