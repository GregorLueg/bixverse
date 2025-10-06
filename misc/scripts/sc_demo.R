# demo for the single cell approach on an mtx file -----------------------------

devtools::load_all()

## object generation -----------------------------------------------------------

# generate the object

object = single_cell_exp(dir_data = tempdir())

# load in the mtx file

# several things are happening here...
# 1.) we scan which genes to include
# 2.) we scan which cells to include
# 3.) we load in the data and save to a binarised CSR file and generate at the
# same time the normalised counts
# 4.) we take the raw counts and normalised counts from the cells and save
# them additionally in a CSC format again to disk.

tictoc::tic()
object = load_mtx(
  object = object,
  sc_mtx_io_param = params_sc_mtx_io(
    path_mtx = path.expand("~/Downloads/ex053/DGE.mtx"),
    path_obs = path.expand("~/Downloads/ex053/cell_metadata.csv"),
    path_var = path.expand("~/Downloads/ex053/all_genes.csv"),
    cells_as_rows = TRUE
  ),
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 100L,
    min_lib_size = 250L,
    min_cells = 10L
  )
)
tictoc::toc()

# show size of the object
pryr::object_size(object)

# seurat takes ~1 minute to do this; scanpy ~2 minutes

list.files(tempdir())

## working with the object -----------------------------------------------------

# get obs table
object[[]]

# get specific columns obs table
object[[c("cell_idx", "cell_id")]]


# get vars
get_sc_var(object)

# demo count retrieval

# gene centric version
tictoc::tic()
counts_gene <- object[, 1:20L, assay = "norm", return_format = "gene"]
tictoc::toc()

counts_gene

class(counts_gene)

tictoc::tic()
counts_gene <- object[, 1:20L, assay = "norm", return_format = "cell"]
tictoc::toc()

class(counts_gene)

object[1:20L, , assay = "raw", return_format = "cell"]

## calculate proportions of some genes of interest -----------------------------

# get gene set proportions
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

G_list <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol"),
  values = get_gene_names(object),
  mart = mart
)
setDT(G_list)

mt_genes <- G_list[external_gene_name %like% "MT-", ensembl_gene_id]
rps_genes <- G_list[external_gene_name %like% "^RPS", ensembl_gene_id]

gs_of_interest <- list(
  "mt_perc" = mt_genes,
  "rb_perc" = rps_genes
)

object <- gene_set_proportions_sc(object, gs_of_interest)

mito_pcs <- object[["mt_perc"]]

hist(unlist(mito_pcs), breaks = 25L)

cells_to_keep <- object[[]][mt_perc <= 0.1, cell_id]

object <- set_cell_to_keep(object, cells_to_keep)

## highly variable genes and pca -----------------------------------------------

# identify HVGs
object <- find_hvg_sc(object = object)

# run PCA
object <- calculate_pca_sc(object, no_pcs = 30L)

## knn graph generation --------------------------------------------------------

object <- find_neighbours_sc(
  object
)

object <- find_clusters_sc(object, res = 1.5)

clusters <- object[["leiden_clustering"]]

samples_original <- unlist(object[["sample"]])

samples_updated <- data.table::tstrsplit(samples_original, "_")
names(samples_updated) <- c("sample_id", "time_point")

object[[names(samples_updated)]] <- samples_updated

samples <- unlist(object[["sample_id"]])

table(unlist(clusters), samples)

demuxlet_data <- fread("~/Downloads/demuxlet_result.best")[,
  c("classification", "classified_sample") := {
    split_result <- strsplit(BEST, "-", fixed = TRUE)
    list(
      sapply(split_result, `[`, 1),
      sapply(split_result, function(x) paste(x[-1], collapse = "-"))
    )
  }
]

filtered_cells <- get_cell_names(object, filtered = TRUE)

doublets <- demuxlet_data[classification == "DBL", BARCODE]

demuxlet_data_red <- demuxlet_data[
  BARCODE %in% get_cell_names(object, filtered = TRUE) & classification == "SNG"
]

head(demuxlet_data_red)

# make it 1 -index
knn_data <- get_knn_mat(object)
rownames(knn_data) <- get_cell_names(object, filtered = TRUE)

knn_edge_list <- rs_knn_mat_to_edge_list(knn_data, TRUE)

label_vec <- data.table(cell_id = rownames(knn_data)) %>%
  merge(
    .,
    demuxlet_data_red[, c("BARCODE", "classified_sample")],
    by.x = "cell_id",
    by.y = "BARCODE",
    all.x = TRUE
  )

predicted_samples <- knn_graph_label_propagation(
  edge_list = knn_edge_list,
  labels = label_vec$classified_sample
)

doublets_bool <- rownames(knn_data) %in% doublets

table(list(
  predictions = predicted_samples$final_labels[
    !doublets_bool &
      is.na(
        label_vec$classified_sample
      )
  ],
  actual = samples[
    !doublets_bool &
      is.na(
        label_vec$classified_sample
      )
  ]
))

table(list(
  predictions = predicted_samples$final_labels[is.na(
    label_vec$classified_sample
  )],
  actual = samples[is.na(label_vec$classified_sample)]
))


f1_score_confusion_mat(
  predicted_samples$final_labels[is.na(
    label_vec$classified_sample
  )],
  samples[is.na(label_vec$classified_sample)]
)

f1_score_confusion_mat(
  predicted_samples$final_labels[
    !doublets_bool &
      is.na(
        label_vec$classified_sample
      )
  ],
  samples[
    !doublets_bool &
      is.na(
        label_vec$classified_sample
      )
  ]
)

f1_score_confusion_mat <- function(clusters_a, clusters_b) {
  # checks
  len_a <- length(clusters_a)
  len_b <- length(clusters_b)
  checkmate::qassert(clusters_a, sprintf("A%i", len_b))
  checkmate::qassert(clusters_b, sprintf("A%i", len_a))

  # function
  cm <- table(clusters_a, clusters_b)

  best_match <- apply(cm, 1, which.max)

  f1_scores <- sapply(seq_len(nrow(cm)), function(i) {
    tp <- cm[i, best_match[i]]
    fp <- sum(cm[, best_match[i]]) - tp
    fn <- sum(cm[i, ]) - tp

    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)

    2 * precision * recall / (precision + recall)
  })

  names(f1_scores) <- rownames(cm)

  return(f1_scores)
}
