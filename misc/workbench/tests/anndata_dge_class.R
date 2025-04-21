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
devtools::check()

object <- bulk_dge_from_h5ad(h5_file) %>%
  change_gene_identifier(alternative_gene_id = "Symbol") %>%
  calculate_pca_bulk_dge()

plot_pca_res(object)

genes_to_filter <- object@variable_info[GeneType == 'protein-coding', Symbol]

object <- calculate_all_dges(
  object = object,
  contrast_column = 'contrast_info',
  filter_column = 'sample_source',
  gene_filter = genes_to_filter
)

get_outputs(object)

# DGE code ----

meta_data <- object@meta_data
dge_list <- object@outputs$dge_list
main_contrast <- 'contrast_info'
co_variates <- NULL
source_col <- 'sample_source'
.verbose <- TRUE

unique(meta_data$sample_source)

## Limma Voom ----

meta_data_red <- meta_data[sample_source == "whole blood"]
dge_list_red <- dge_list[genes_to_filter, meta_data_red$sample_id]

main_contrast <- "contrast_info"

?run_limma_voom

dge_results <- run_limma_voom(
  meta_info = meta_data_red,
  main_contrast = main_contrast,
  dge_list = dge_list_red
)

## Effect sizes ----

effect_size_results <- hedges_g_dge_list(
  meta_info = meta_data_red,
  main_contrast = main_contrast,
  dge_list = dge_list_red
)
