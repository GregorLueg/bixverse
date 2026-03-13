vars <- get_sc_var(single_cell_obj)

tfs <- fread(
  "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt",
  header = FALSE
) %>%
  `colnames<-`("tf")

tf_indices <- which(vars$column1 %in% tfs$tf) - 1L

object = single_cell_obj

genes_to_include <- rs_scenic_gene_filter(
  f_path_genes = get_rust_count_gene_f_path(object),
  cell_indices = get_cells_to_keep(object),
  scenic_params = params_scenic(min_counts = 15L),
  verbose = TRUE
)

includes_genes_ids <- vars[genes_to_include + 1, column1]

tf_indices_red <- intersect(tf_indices, genes_to_include)

scenic_results <- rs_scenic_grn(
  f_path_genes = get_rust_count_gene_f_path(object),
  cell_indices = get_cells_to_keep(object),
  gene_indices = genes_to_include,
  tf_indices = as.integer(tf_indices_red),
  scenic_params = params_scenic(
    gene_batch_size = 64L,
    learner_type = "randomforest"
  ),
  seed = 123L,
  verbose = TRUE
)

top_x_per_column <- function(mat, x, gene_ids) {
  apply(
    mat,
    2,
    function(col) gene_ids[order(col, decreasing = TRUE)[1:x]],
    simplify = FALSE
  )
}

results <- top_x_per_column(
  scenic_results,
  x = 50L,
  gene_ids = includes_genes_ids
)

names(results) <- vars$column1[tf_indices_red + 1]

paths <- download_cistarget_hg38()
rankings <- read_motif_ranking(paths$rankings)
annotations <- read_motif_annotation_file(paths$motif_annotations)

cis_target_results <- run_cistarget(
  gs_list = results,
  rankings = rankings,
  annot_data = annotations
)

cis_target_results[,
  match := purrr::map2_lgl(
    gs_name,
    TF_highConf,
    \(tf, associated_tf) {
      if (is.na(associated_tf)) {
        return(FALSE)
      }
      tf %in% strsplit(associated_tf, ";")[[1]]
    }
  )
]

cis_target_results[(match)]$gs_name
