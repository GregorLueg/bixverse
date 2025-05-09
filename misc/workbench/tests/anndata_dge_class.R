library("anndata")
library(data.table)
library(magrittr)

library(reticulate)

reticulate::py_list_packages()

reticulate::py_install('anndata')

ad <- read_h5ad("~/Desktop/geo_data/GSE154773_anndata.h5ad")

rownames(ad$obs) == h5_obj$get_obs_table()$sample_id

all(rownames(ad$var) == h5_obj$get_var_info()$var_id)

all()

counts <- t(ad$X)

devtools::load_all()
devtools::document()

obs_original = ad$obs

h5_file <- "~/Desktop/geo_data/GSE154773_anndata.h5ad"

h5_obj <- anndata_parser$new(h5_file)

all(rownames(test_counts) == h5_obj$get_var_info()$var_id)

test_counts <- h5_obj$get_raw_counts()

test_counts["3553", ]

c(meta_data, var_info, counts) %<-% h5_obj$get_key_data()

counts["3553", ]

all(rownames(counts) == rownames(test_counts))

object <- bulk_dge(
  raw_counts = counts,
  meta_data = meta_data,
  variable_info = var_info
)

all(rownames(object@raw_counts) == var_info$var_id)

IL1B_original <- object@raw_counts["3553", ]

object <- change_gene_identifier(object, "Symbol")

object@raw_counts["IL1B", ]

object@outputs$dge_list["IL1B", ]

test_counts[1:5, 1:5]

if (.verbose) message("Loading data from the h5ad object")
c(meta_data, var_info, counts)
h5_contentc(meta_data, var_info, counts) %<-% h5_obj$get_key_data()


h5_obj$get_var_info()

bulk_dge(raw_counts = counts, meta_data = meta_data, variable_info = var_info)

devtools::load_all()
devtools::document()
devtools::check()

object <- bulk_dge_from_h5ad(h5_file)
object <- change_gene_identifier(object, "Symbol")
object <- calculate_pca_bulk_dge(object)

object <- calculate_all_dges(
  object = object,
  contrast_column = 'contrast_info',
  filter_column = 'sample_source',
  gene_filter = genes_to_filter
)


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

# Test crazy bug ----

mat <- matrix(data = rnorm(10 * 10), ncol = 10, nrow = 10)
mat2 <- matrix(data = rnorm(10 * 10), ncol = 10, nrow = 10)
rownames(mat) <- rownames(mat2) <- sprintf("test_%i", 10:1)

?copy

actual_var <- copy(rownames(mat))

weird_dt <- list(
  value = actual_var
) %>%
  setDT() %>%
  setorder(value)
