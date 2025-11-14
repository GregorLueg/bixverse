# test seacell -----------------------------------------------------------------

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")

rextendr::document()

devtools::load_all()

sessionInfo()

data_dir <- path.expand("~/Desktop/seacell_test/")

sc_object <- single_cell_exp(dir_data = tempdir())

h5_path <- path.expand(
  "~/repos/other/Hotspot/notebooks/5k_h5ad.h5ad"
)

sc_object <- load_h5ad(
  sc_object,
  h5_path = h5_path,
  streaming = FALSE
)

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = 2000L,
  .verbose = FALSE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 30L,
  randomised_svd = TRUE,
  .verbose = FALSE
)

# sea_cells <- get_seacells_sc(
#   sc_object,
#   seacell_params = params_sc_seacells(
#     n_sea_cells = 90L,
#     knn = list(k = 15L),
#     convergence_epsilon = 1e-5,
#     max_fw_iters = 50L,
#     pruning = TRUE
#   )
# )

object = sc_object

hotspot_params = params_sc_hotspot(
  model = "danb",
  knn = list(ann_dist = "cosine", k = 30L),
  normalise = TRUE
)

tictoc::tic()
hotspot_auto_cor <- rs_hotspot_autocor(
  f_path_genes = get_rust_count_gene_f_path(object),
  f_path_cells = get_rust_count_cell_f_path(object),
  embd = get_pca_factors(object),
  hotspot_params = hotspot_params,
  cells_to_keep = get_cells_to_keep(object),
  genes_to_use = get_gene_indices(
    object,
    gene_ids = get_gene_names(object),
    rust_index = TRUE
  ),
  streaming = FALSE,
  verbose = TRUE,
  seed = 42L
)
tictoc::toc()

hotspot_auto_cor_dt <- as.data.table(hotspot_auto_cor)[,
  gene_name := get_gene_names(object)
]

setorder(hotspot_auto_cor_dt, -z_score)

hotspot_python <- fread(
  "~/repos/other/Hotspot/notebooks/hotspot_python_results.csv"
)

combined = merge(
  hotspot_python[, c("Gene", "C", "Z")],
  hotspot_auto_cor_dt[, c("gaerys_c", "gene_name", "z_score")],
  by.x = "Gene",
  by.y = "gene_name"
)

plot(combined$C, combined$gaerys_c)

plot(combined$Z, combined$z_score)

tictoc::tic()
hotspot_gene_cor <- rs_hotspot_gene_cor(
  f_path_genes = get_rust_count_gene_f_path(object),
  f_path_cells = get_rust_count_cell_f_path(object),
  embd = get_pca_factors(object),
  hotspot_params = hotspot_params,
  cells_to_keep = get_cells_to_keep(object),
  genes_to_use = as.integer(hotspot_auto_cor_dt$gene_idx[1:2500]),
  streaming = TRUE,
  verbose = TRUE,
  seed = 42L
)
tictoc::toc()

heatmap(hotspot_gene_cor$z)

hotspot_gene_cor$z[1:10, 1:10]

hist(hotspot_auto_cor$gaerys_c)

hist(hotspot_auto_cor$z_score)

table(hotspot_auto_cor$fdr <= 0.05)

seacell_result <- fread("/Users/gregor/Desktop/seacell_results.csv")

class(object)

plot(
  x = seq_len(length(sea_cells@other_data$rss)),
  y = sea_cells@other_data$rss
)


object = sea_cells
original_cell_type = seacell_result$celltype


sea_cells <- calc_meta_cell_purity(sea_cells, seacell_result$celltype)


internal_results <- data.table(
  cell_id = sc_object[["index"]],
  assignment = sea_cells@original_assignment$assignments
)[, cell_type := seacell_result$celltype]


table(internal_results$assignment)

mclust::adjustedRandIndex(internal_results$assignment, seacell_result$SEACell)

sea_cells@obs_table$no_originating_cells

# > sea_cells@obs_table$no_originating_cells
#  [1]  56  48  47  78  69  59 127  90  66 210  34  49  35  50  66  93  76  37  40  36 180  43  65  76  54 246  71  34  81  66  89  58  49  47  41  48  34
# [38]  25  93  31  72  64  46  47  82 166 102 234  60  40  33  34 105  49  65 123  34 169 119 164  19  34  71  32  23 150  76  93  53  49  83  20 112  87
# [75] 124  71  72  48 128 210  19  32  82 155 178  22  58  62  42  71

# > sea_cells@obs_table$no_originating_cells
#  [1]  57  51  49  79  69  56 119  87  68 178  33  49  37  53  66  93  77  38  43  41 170  46  63  69  55 247  76  31  74  67  87  56  45  50  48  45  35
# [38]  24 101  29  69  70  46  42  80 166 105 263  60  40  35  34 100  50  59 141  35 177 120 158  20  33  71  35  23 150  79 100  47  42  78  21 118  81
# [75] 127  80  69  52 128 202  19  30  82 161 169  18  56  62  42  75

rss_vals_wo_pruning <- sea_cells@other_data$rss
