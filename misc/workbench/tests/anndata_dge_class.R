library("anndata")
library(data.table)
library(magrittr)

library(reticulate)

reticulate::py_list_packages()

reticulate::py_install('anndata')

ad <- read_h5ad("~/Desktop/geo_data/GSE154773_anndata.h5ad")

all(rownames(ad$var) == h5_obj$get_var_info()$var_id)

counts <- t(ad$X)

devtools::load_all()
devtools::document()

obs_original = ad$obs

h5_file <- "~/Desktop/geo_data/GSE154773_anndata.h5ad"

h5_path <- h5_file

h5_obj <- anndata_parser$new(h5_file)


if (.verbose) message("Loading data from the h5ad object")
c(meta_data, var_info, counts)
h5_contentc(meta_data, var_info, counts) %<-% h5_obj$get_key_data()


h5_obj$get_var_info()

bulk_dge(raw_counts = counts, meta_data = meta_data, variable_info = var_info)

devtools::load_all()
devtools::document()


object <- bulk_dge_from_h5ad(h5_file) %>%
  change_gene_identifier(alternative_gene_id = "Symbol") %>%
  calculate_pca_bulk_dge()

plot_pca_res(object)

# DGE code ----

genes_to_filter <- object@variable_info[GeneType == 'protein-coding', Symbol]
meta_data <- object@meta_data
dge_list <- object@outputs$dge_list
main_contrast <- 'contrast_info'
co_variates <- NULL
source_col <- 'sample_source'
.verbose <- TRUE

unique(meta_data$sample_source)

## Limma Voom ----

meta_data_red <- meta_data[sample_source == "skin"]
dge_list_red <- dge_list[genes_to_filter, meta_data_red$sample_id]

main_contrast <- "contrast_info"

dge_results <- run_limma_voom(
  meta_info = meta_data_red,
  main_contrast = main_contrast,
  dge_list = dge_list_red
)

## Effect sizes ----

groups <- as.character(unique(meta_data_red[['contrast_info']]))

combinations_to_test <- combn(
  x = groups,
  m = 2,
  FUN = function(x) {
    c(x[[1]], x[[2]])
  },
  simplify = FALSE
)


combination <- combinations_to_test[[1]]

res <- purrr::map(combinations_to_test, \(combination) {
  grpA <- meta_data_red[
    eval(parse(text = paste0(main_contrast, " == '", combination[[1]], "'"))),
    sample_id
  ]
  grpB <- meta_data_red[
    eval(parse(text = paste0(main_contrast, " == '", combination[[2]], "'"))),
    sample_id
  ]

  to_keep <- suppressWarnings(edgeR::filterByExpr(dge_list_red, design = NULL))

  voom_obj <- limma::voom(counts = dge_list_red[to_keep, ])

  mat_a <- t(voom_obj$E[, grpA])
  mat_b <- t(voom_obj$E[, grpB])

  hedges_g_effect <- calculate_effect_size(mat_a = mat_a, mat_b = mat_b) %>%
    data.table::setDT() %>%
    .[, `:=`(
      gene_id = colnames(mat_a),
      combination = paste(combination[[1]], combination[[2]], sep = "_vs_")
    )]

  hedges_g_effect
})


sum(to_keep)

model_formula

rm(model_formula)
