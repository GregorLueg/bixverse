# atac on the multi-modal class -------------------------------------------------

# port of bixverse-rs python/tests/test_atac.py: numbers against Signac, MACS3
# and AMULET on the 500 PBMC multiome cells (chr22)

source("helper_atac.R", local = TRUE)
if (!ATAC_OK) {
  exit_file("ATAC fixtures not found; set BIXVERSE_ATAC_FIXTURES")
}

## store and qc ----------------------------------------------------------------

mm <- atac_mm_experiment()
signac <- atac_tsv("signac/cell_metrics.tsv")

obs <- get_sc_obs(mm, cols = c("cell_id", "atac_lib_size"))
expect_equal(obs$cell_id, signac$barcode, info = "store in barcode order")
expect_equal(
  as.integer(obs$atac_lib_size),
  signac$nCount_peaks,
  info = "fragment counts match Signac nCount_peaks"
)
var <- get_sc_var(mm, modality = "atac")
expect_equal(nrow(var), 1883L, info = "1883 peaks")
expect_true(startsWith(var$feature_id[1], "chr22:"), info = "region ids")

mm <- fragment_qc_sc(mm, ATAC_FRAGS, ATAC_GTF)
obs <- get_sc_obs(mm)
expect_equal(
  as.integer(obs$atac_n_fragments),
  signac$frequency_count,
  info = "fragments per cell match Signac"
)
expect_equal(
  obs$atac_tss_enrichment,
  signac$tss_enrichment_fast,
  tolerance = 1e-9,
  info = "TSS enrichment matches Signac"
)
ok <- !is.na(signac$frip)
expect_equal(
  obs$atac_frip[ok],
  signac$frip[ok],
  tolerance = 1e-12,
  info = "FRiP matches Signac"
)

## tf-idf, top features and lsi ------------------------------------------------

expect_error(
  calculate_lsi_sc(mm, lsi_params = params_atac_lsi(n_comp = 5L)),
  pattern = "run_tfidf",
  info = "LSI needs TF-IDF first"
)

mm <- run_tfidf_sc(mm)
mm <- find_top_features_sc(mm, min_cutoff = "q5")
ref_top <- atac_tsv("signac/top_features.tsv")
expect_equal(
  which(get_sc_var(mm, modality = "atac")$top_feature) - 1L,
  which(ref_top$kept == 1) - 1L,
  info = "top features match Signac"
)

mm <- calculate_lsi_sc(
  mm,
  lsi_params = params_atac_lsi(
    n_comp = 10L,
    drop_first = FALSE,
    solver = "exact"
  ),
  .verbose = FALSE
)
emb <- get_embedding(mm, "lsi", modality = "atac")
ref_emb <- atac_tsv("signac/lsi_embeddings.tsv.gz")
expect_equal(dim(emb), c(500L, 10L), info = "all 10 components kept")
cors <- purrr::map_dbl(seq_len(10L), \(j) {
  abs(stats::cor(emb[, j], ref_emb[[sprintf("LSI_%i", j)]]))
})
expect_true(all(cors > 0.999), info = "LSI matches Signac up to sign")
ref_depth <- atac_tsv("signac/lsi_components.tsv")$depth_cor
expect_equal(
  abs(get_lsi_depth_cor(mm)$depth_cor),
  abs(ref_depth),
  tolerance = 1e-4,
  info = "depth correlation matches Signac"
)

mm <- calculate_lsi_sc(
  mm,
  lsi_params = params_atac_lsi(n_comp = 10L),
  .verbose = FALSE
)
expect_equal(
  dim(get_embedding(mm, "lsi", modality = "atac")),
  c(500L, 9L),
  info = "drop_first leaves 9 components"
)
expect_equal(
  get_lsi_depth_cor(mm)$dropped[1:2],
  c(TRUE, FALSE),
  info = "component 1 is flagged as dropped"
)

mm <- find_hvg_sc(mm, hvg_no = 40L, .verbose = FALSE)
mm <- calculate_pca_sc(mm, no_pcs = 10L, .verbose = FALSE)
mm <- generate_wnn_graph_sc(
  mm,
  modality_2 = "atac",
  embd_to_use_2 = "lsi",
  .verbose = FALSE
)
expect_equal(
  nrow(get_knn_mat(mm, modality = "wnn")),
  500L,
  info = "WNN over RNA and the ATAC LSI"
)

## errors ----------------------------------------------------------------------

bare <- atac_mm_experiment(peaks = FALSE)
expect_error(run_tfidf_sc(bare), pattern = "No ATAC", info = "needs a store")
expect_error(
  find_markers_sc(
    bare,
    ATAC_BARCODES[1:10],
    ATAC_BARCODES[11:20],
    modality = "atac",
    .verbose = FALSE
  ),
  pattern = "No ATAC",
  info = "differential accessibility needs a store"
)

## tiles -----------------------------------------------------------------------

tiles <- add_atac_fragments_sc(
  bare,
  ATAC_FRAGS,
  tile_width = 5000L,
  chrom_sizes = c(chr22 = 50818468L),
  mode = "fragment",
  .verbose = FALSE
)
expect_equal(tiles@atac_info$kind, "tile", info = "tile store")
expect_equal(
  nrow(get_sc_var(tiles, modality = "atac")),
  10164L,
  info = "5 kb tiles over chr22"
)
expect_error(
  add_atac_fragments_sc(tiles, ATAC_FRAGS, overwrite = TRUE),
  pattern = "exactly one",
  info = "peaks or tiles, not neither"
)
expect_error(
  find_top_features_sc(tiles, min_cutoff = "x5"),
  pattern = "min_cutoff",
  info = "malformed cutoff"
)

## peak calling ----------------------------------------------------------------

pc <- atac_mm_experiment(peaks = FALSE)
get_sc_duckdb(pc)$add_data_obs(data.table::data.table(
  macs_group = sprintf("g%i", (seq_len(500L) - 1L) %% 3L)
))
per_group <- call_peaks_sc(
  pc,
  ATAC_FRAGS,
  group_by = "macs_group",
  combine = FALSE
)
for (g in c("g0", "g1", "g2")) {
  ref <- atac_tsv(sprintf("macs3/%s_insertion_dupall_g2.7e9.tsv.gz", g))
  ours <- per_group[as.character(per_group$group) == g]
  for (col in c("chrom", "start", "end", "summit")) {
    expect_equal(
      as.vector(ours[[col]]),
      as.vector(ref[[col]]),
      info = sprintf("MACS3 %s %s", g, col)
    )
  }
  expect_equal(
    ours$qscore,
    ref$qscore,
    tolerance = 1e-5,
    info = sprintf("MACS3 %s qscore", g)
  )
}

signac_style <- call_peaks_sc(
  pc,
  ATAC_FRAGS,
  combine = FALSE,
  peak_params = params_atac_call_peaks(
    tag_source = "fragment_start",
    keep_dup = 1L,
    genome_size = 1e7
  )
)
ref <- atac_tsv("macs3/all_start_dup1_g1e7.tsv.gz")
expect_equal(
  signac_style$start,
  as.integer(ref$start),
  info = "Signac's MACS call"
)
expect_equal(
  unique(as.character(signac_style$group)),
  "all",
  info = "pooled group"
)

merged <- call_peaks_sc(pc, ATAC_FRAGS, group_by = "macs_group")
expect_equal(
  names(merged),
  c("chrom", "start", "end", "peak_called_in"),
  info = "combined peaks"
)
expect_true(
  nrow(merged) > 0 && nrow(merged) < nrow(per_group),
  info = "merging shrinks the peak set"
)
expect_true(
  all(unlist(strsplit(merged$peak_called_in, ",")) %in% c("g0", "g1", "g2")),
  info = "called-in groups"
)
pc <- add_atac_fragments_sc(pc, ATAC_FRAGS, peaks = merged, .verbose = FALSE)
expect_equal(
  nrow(get_sc_var(pc, modality = "atac")),
  nrow(merged),
  info = "called peaks feed counting"
)

## amulet ----------------------------------------------------------------------

am <- amulet_sc(
  bare,
  ATAC_FRAGS,
  restrict_beds = file.path(ATAC_FIX, "blacklist.bed")
)
ref <- atac_tsv("amulet/overlap2_MultipletProbabilities.txt.gz")
obs <- get_sc_obs(
  am,
  cols = c("cell_id", "amulet_p", "amulet_q", "amulet_multiplet")
)
joined <- merge(obs, ref, by = "cell_id")
expect_equal(nrow(joined), 500L, info = "every cell tested")
expect_equal(
  joined$amulet_p,
  joined$`p-value`,
  tolerance = 1e-9,
  info = "AMULET p-values"
)
expect_equal(
  joined$amulet_q,
  joined$`q-value`,
  tolerance = 1e-9,
  info = "AMULET q-values"
)
expect_true(
  sum(joined$amulet_multiplet) == sum(joined$`q-value` < 0.01) &&
    sum(joined$amulet_multiplet) > 0,
  info = "multiplet calls"
)

## gene activity ---------------------------------------------------------------

ga <- gene_activity_sc(bare, ATAC_FRAGS, ATAC_GTF, .verbose = FALSE)
expect_equal(
  rs_store_feature_kind(bixverse:::get_rust_count_cell_f_path(ga)),
  "gene_activity",
  info = "gene activity store kind"
)
expect_equal(
  rs_store_feature_kind(bixverse:::get_rust_count_cell_f_path(bare)),
  "gene",
  info = "the RNA stays genes"
)
expect_equal(
  get_cell_names(ga),
  get_cell_names(bare),
  info = "same cells, same order"
)

names_ref <- atac_tsv("signac/gene_activity/features.tsv.gz")$name
ref <- atac_tsv("signac/gene_activity/matrix.tsv.gz")
var <- get_sc_var(ga)
expect_true(
  startsWith(var$gene_id[1], "ENSG") &&
    !grepl(".", var$gene_id[1], fixed = TRUE),
  info = "unversioned Ensembl ids"
)
expect_equal(sort(var$gene_name), sort(names_ref), info = "Signac's genes")
counts <- as(get_sc_counts(ga, assay = "raw"), "TsparseMatrix")
got <- sort(sprintf(
  "%s|%i|%i",
  var$gene_name[counts@j + 1L],
  counts@i,
  as.integer(counts@x)
))
want <- sort(sprintf(
  "%s|%i|%i",
  names_ref[ref$feature_idx + 1L],
  ref$cell_idx,
  as.integer(ref$value)
))
expect_equal(got, want, info = "gene activity counts match Signac")
expect_equal(
  get_sc_obs(ga)$lib_size,
  as.numeric(Matrix::rowSums(get_sc_counts(ga, assay = "raw"))),
  info = "lib_size is the gene activity's"
)

expect_error(
  scrublet_sc(ga),
  pattern = "needs RNA counts",
  info = "scrublet refuses gene activity"
)
expect_warning(
  bixverse:::.check_rna_counts(ga, "aucell", refuse = FALSE),
  pattern = "holds gene activity",
  info = "AUCell warns on gene activity"
)
expect_silent(bixverse:::.check_rna_counts(bare, "scrublet", refuse = TRUE))

## frip in insertion mode ------------------------------------------------------

ins <- add_atac_fragments_sc(
  bare,
  ATAC_FRAGS,
  peaks = ATAC_PEAKS,
  overwrite = TRUE,
  .verbose = FALSE
)
expect_equal(ins@atac_info$count_mode, "insertion", info = "default mode")
ins <- fragment_qc_sc(ins, ATAC_FRAGS, ATAC_GTF)
obs <- get_sc_obs(ins)
ok <- obs$atac_n_fragments > 0
expect_equal(
  obs$atac_frip[ok],
  obs$atac_lib_size[ok] / (2 * obs$atac_n_fragments[ok]),
  tolerance = 1e-12,
  info = "insertion FRiP"
)
expect_true(max(obs$atac_frip[ok]) <= 1, info = "insertion FRiP <= 1")

## differential accessibility --------------------------------------------------

da <- atac_mm_experiment()
m <- atac_tsv("signac/peaks/matrix.tsv.gz")
dense <- matrix(0, 1883L, 500L)
dense[cbind(m$feature_idx + 1L, m$cell_idx + 1L)] <- m$value
groups <- list(a = ATAC_BARCODES[1:250], b = ATAC_BARCODES[251:500])

pb <- get_pseudobulked_sc(da, groups, modality = "atac", .verbose = FALSE)
expect_equal(
  unname(pb),
  rbind(rowSums(dense[, 1:250]), rowSums(dense[, 251:500])),
  info = "ATAC pseudobulk"
)
expect_equal(
  colnames(pb),
  get_sc_var(da, modality = "atac")$feature_id,
  info = "pseudobulk keyed by feature id"
)

res <- find_markers_sc(
  da,
  groups$a,
  groups$b,
  modality = "atac",
  .verbose = FALSE
)
expect_true(startsWith(res$feature_id[1], "chr22:"), info = "feature ids")
lr <- find_markers_sc(
  da,
  groups$a,
  groups$b,
  method = "lr",
  modality = "atac",
  .verbose = FALSE
)
expect_equal(res, lr, info = "ATAC defaults to the LR test")
wil <- find_markers_sc(
  da,
  groups$a,
  groups$b,
  method = "wilcox",
  modality = "atac",
  .verbose = FALSE
)
expect_false(
  isTRUE(all.equal(wil$p_values, res$p_values)),
  info = "Wilcoxon differs from LR"
)
arms <- find_all_markers_sc(
  da,
  "cell_grp",
  modality = "atac",
  .verbose = FALSE
)
expect_true(startsWith(arms$feature_id[1], "chr22:"), info = "all markers")
spec <- find_specific_markers_sc(
  da,
  "cell_grp",
  modality = "atac",
  .verbose = FALSE
)
expect_equal(spec$params$method, "lr", info = "specific markers default to LR")

## harmony on the lsi ----------------------------------------------------------

hm <- atac_mm_experiment()
hm <- run_tfidf_sc(hm)
hm <- find_top_features_sc(hm, min_cutoff = "q5")
hm <- calculate_lsi_sc(
  hm,
  lsi_params = params_atac_lsi(n_comp = 10L),
  .verbose = FALSE
)
get_sc_duckdb(hm)$add_data_obs(data.table::data.table(
  atac_batch = ifelse(seq_len(500L) %% 2L == 1L, "a", "b")
))
expect_warning(
  harmony_sc(hm, "atac_batch", modality = "atac", .verbose = FALSE),
  pattern = "'pca'",
  info = "no ATAC PCA"
)
hm <- harmony_sc(
  hm,
  "atac_batch",
  embd_to_use = "lsi",
  modality = "atac",
  .verbose = FALSE
)
lsi <- get_embedding(hm, "lsi", modality = "atac")
corrected <- get_embedding(hm, "harmony", modality = "atac")
expect_equal(dim(corrected), dim(lsi), info = "harmony keeps the shape")
expect_false(isTRUE(all.equal(unname(corrected), unname(lsi))), info = "moved")
expect_false(
  "harmony" %in% get_available_embeddings(hm, modality = "rna"),
  info = "written to the ATAC cache only"
)
hm <- find_neighbours_sc(
  hm,
  embd_to_use = "harmony",
  modality = "atac",
  .verbose = FALSE
)
expect_equal(
  nrow(get_knn_mat(hm, modality = "atac")),
  500L,
  info = "kNN on the corrected LSI"
)
