# fix hotspot ------------------------------------------------------------------

devtools::load_all()

dir_data <- path.expand("~/Downloads/single_cell_data/")

sc_object <- single_cell_exp(
  dir_data = dir_data
)

sc_object <- load_existing(sc_object)

dim(get_knn_mat(sc_object))

dim(get_pca_factors(sc_object))

hotspot_autocor_danb_res <- hotspot_autocor_sc(
  hotspot_params = params_sc_hotspot(),
  object = sc_object,
  embd_to_use = "pca",
  streaming = TRUE,
  .verbose = TRUE
)

get_sc_var(sc_object)

weird_genes <- get_gene_indices(
  sc_object,
  gene_ids = hotspot_autocor_danb_res$gene_id[is.na(
    hotspot_autocor_danb_res$pval
  )],
  rust_index = FALSE
)

Matrix::colSums(sc_object[, weird_genes, return_format = "gene"])
