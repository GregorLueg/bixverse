# test seacell -----------------------------------------------------------------

rextendr::document()

data_dir <- path.expand("~/Desktop/seacell_test/")

sc_object <- single_cell_exp(dir_data = data_dir)

h5_path <- path.expand(
  "~/repos/other/SEACells/notebooks/data/cd34_multiome_rna.h5ad"
)

sc_object <- load_h5ad(
  sc_object,
  h5_path = path.expand(
    "~/repos/other/SEACells/notebooks/data/cd34_multiome_rna.h5ad"
  ),
  streaming = FALSE
)

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = 1500L,
  .verbose = FALSE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 50L,
  randomised_svd = TRUE,
  .verbose = FALSE
)

sea_cells <- get_seacells_sc(
  sc_object,
  seacell_params = params_sc_seacells(
    n_sea_cells = 90L,
    k = 15L,
    convergence_epsilon = 1e-5,
    max_fw_iters = 50L,
    pruning = TRUE
  )
)

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


# test on other data -----------------------------------------------------------

devtools::load_all()

dir_data <- path.expand(
  "~/Downloads/splitpipe/Mosaic_POC_Combined/all-sample/DGE_unfiltered/"
)

list.files(dir_data)

sc_object <- single_cell_exp(dir_data = tempdir())

sc_object <- load_mtx(
  sc_object,
  sc_mtx_io_param = params_sc_mtx_io(
    path_mtx = file.path(dir_data, "count_matrix.mtx"),
    path_obs = file.path(dir_data, "cell_metadata.csv"),
    path_var = file.path(dir_data, "all_genes.csv"),
    cells_as_rows = TRUE,
    has_hdr = TRUE
  ),
  sc_qc_param = params_sc_min_quality()
)

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = 2000L,
  .verbose = TRUE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 30L,
  randomised_svd = TRUE,
  .verbose = TRUE
)

sea_cells <- get_seacells_sc(
  sc_object,
  seacell_params = params_sc_seacells(
    n_sea_cells = 1500L,
    k = 15L,
    convergence_epsilon = 1e-5,
    max_fw_iters = 50L,
    pruning = TRUE
  )
)
