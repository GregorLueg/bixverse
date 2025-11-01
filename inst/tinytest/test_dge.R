# dge --------------------------------------------------------------------------

library(magrittr)

## data ------------------------------------------------------------------------

### test data ------------------------------------------------------------------

test_data <- qs2::qs_read(
  "./synthetic_data/dge_test_data.qs"
)

### expected data --------------------------------------------------------------

expected_pca_pvals <- c(3.494401e-05, 8.877927e-01)
expected_pc1 <- c(
  37.19602,
  -20.14648,
  22.61833,
  -37.88058,
  30.08617,
  -23.42578,
  27.39613,
  -35.84382
)
expected_pc2 <- c(
  -16.95906,
  -13.90214,
  -11.77074,
  -11.08694,
  40.49074,
  46.80947,
  -17.84254,
  -15.73879
)

expected_limma_res <- data.table::fread("./test_data/dge_limma_res.gz")
data.table::setorder(expected_limma_res, gene_id)

expected_effect_sizes <- data.table::fread("./test_data/dge_effect_sizes.gz")
data.table::setorder(expected_effect_sizes, gene_id)

## test the class --------------------------------------------------------------

### class generation -----------------------------------------------------------

bad_metadata <- data.table::copy(test_data$meta_data) %>%
  data.table::setnames(old = "sample_id", new = "doesn't work")

bad_counts <- test_data$counts
rownames(bad_counts) <- NULL

expect_error(
  current = bulk_dge(
    raw_counts = test_data$counts,
    meta_data = bad_metadata
  ),
  info = paste("DGE class with bad metadata.")
)

expect_error(
  current = bulk_dge(
    raw_counts = bad_counts,
    meta_data = test_data$meta_data
  ),
  info = paste("DGE class with bad counts.")
)

dge_class <- bulk_dge(
  raw_counts = test_data$counts,
  meta_data = test_data$meta_data
)

### qc -------------------------------------------------------------------------

# remove lowly expressed genes and outlier samples
dge_class <- qc_bulk_dge(
  dge_class,
  group_col = "dex",
  min_prop = 0.2,
  min_count = 10,
  .verbose = FALSE
)

qc_plot_1 <- get_dge_qc_plot(dge_class, plot_choice = 1L)
qc_plot_2 <- get_dge_qc_plot(dge_class, plot_choice = 2L)

expect_true(
  current = "ggplot" %in% class(qc_plot_1),
  info = "DGE pre-processing: 1st QC plot"
)
expect_true(
  current = "ggplot" %in% class(qc_plot_2),
  info = "DGE pre-processing: 2nd QC plot"
)

expect_warning(
  current = get_dge_qc_plot(dge_class, plot_choice = 3L),
  info = "Non existing plot 3 warning"
)
expect_warning(
  current = get_dge_qc_plot(dge_class, plot_choice = 4L),
  info = "Non existing plot 4 warning"
)

### normalisations -------------------------------------------------------------

dge_class <- normalise_bulk_dge(dge_class, group_col = "dex", .verbose = FALSE)

qc_plot_3 <- get_dge_qc_plot(dge_class, plot_choice = 3L)
qc_plot_4 <- get_dge_qc_plot(dge_class, plot_choice = 4L)


expect_true(
  current = "ggplot" %in% class(qc_plot_3),
  info = "DGE pre-processing: 3rd QC plot"
)
expect_true(
  current = "ggplot" %in% class(qc_plot_4),
  info = "DGE pre-processing: 4th QC plot"
)
expect_true(
  current = all(dim(get_dge_list(dge_class)) == c(15926, 8)),
  info = "DGE pre-processing - filtered lowly expressed genes."
)
expect_true(
  current = class(get_dge_list(dge_class)) == "DGEList",
  info = "DGE class - get the DGEList object"
)

### pca ------------------------------------------------------------------------

expect_warning(
  current = get_dge_qc_plot(dge_class, plot_choice = 5L),
  info = "Non existing plot 5 warning"
)

dge_class <- calculate_pca_bulk_dge(dge_class)

expect_equal(
  current = dge_class@outputs$pca$PC_1,
  target = expected_pc1,
  tolerance = 1e-7,
  info = "DGE class - PC1 values"
)

expect_equal(
  current = dge_class@outputs$pca$PC_2,
  target = expected_pc2,
  tolerance = 1e-6,
  info = "DGE class - PC2 values"
)

qc_plot_5 <- get_dge_qc_plot(dge_class, plot_choice = 5L)

expect_true(
  current = "ggplot" %in% class(qc_plot_5),
  info = "DGE pre-processing: 5th QC plot"
)

pca_plot_v2 <- plot_pca_res(dge_class, cols_to_plot = c("dex", "cell"))

expect_true(
  current = "ggplot" %in% class(pca_plot_v2),
  info = "DGE pre-processing: 5th QC plot"
)

### dge calculations -----------------------------------------------------------

#### limma ---------------------------------------------------------------------

expect_warning(
  current = get_dge_limma_voom(dge_class),
  info = "DGE class - warning of Limma Voom getter"
)

dge_class <- calculate_dge_limma(
  dge_class,
  contrast_column = "dex",
  .verbose = FALSE
)

limma_res <- get_dge_limma_voom(dge_class)
data.table::setorder(limma_res, gene_id)

expect_equal(
  current = limma_res,
  target = expected_limma_res,
  info = "DGE class - Limma results"
)

#### hedge's g -----------------------------------------------------------------

expect_warning(
  current = get_dge_effect_sizes(dge_class),
  info = "DGE class - warning of Hedge's getter"
)

dge_class <- calculate_dge_hedges(
  dge_class,
  contrast_column = "dex",
  .verbose = FALSE
)

effect_sizes <- get_dge_effect_sizes(dge_class)
data.table::setorder(effect_sizes, gene_id)

expect_equal(
  current = effect_sizes,
  target = expected_effect_sizes,
  info = "DGE class - Effect size results"
)

## normalisation methods -------------------------------------------------------

### gene length data -----------------------------------------------------------

gene_lengths <- c(
  1745.5,
  873.5,
  1074.5,
  2870.0,
  3086.5,
  2383.0,
  3136.0,
  1316.0,
  1800.0,
  3759.0
)
names(gene_lengths) <- rownames(test_data$counts)[1:10]

dge_class <- bulk_dge(
  raw_counts = test_data$counts[1:10, ],
  meta_data = test_data$meta_data
)

dge_class <- qc_bulk_dge(dge_class, group_col = "dex", .verbose = FALSE)

expect_error(
  current = normalise_bulk_dge(
    dge_class,
    group_col = "dex",
    calc_tpm = TRUE,
    .verbose = FALSE
  ),
  info = paste("dge class - error without gene lengths")
)

expect_warning(
  current = get_tpm_counts(dge_class),
  info = paste("dge class - no TPM counts")
)

expect_warning(
  current = get_fpkm_counts(dge_class),
  info = paste("dge class - no FPKM counts")
)

dge_class <- normalise_bulk_dge(
  dge_class,
  group_col = "dex",
  calc_tpm = TRUE,
  calc_fpkm = TRUE,
  gene_lengths = gene_lengths,
  .verbose = FALSE
)

expect_true(
  current = "matrix" %in% class(get_tpm_counts(dge_class)),
  info = paste("dge class - returns TPM counts")
)

expect_true(
  current = "matrix" %in% class(get_fpkm_counts(dge_class)),
  info = paste("dge class - returns FPKM counts")
)
