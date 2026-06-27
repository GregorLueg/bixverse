## Libraries
library(bixverse)
library(bixverse.plots)
library(ggplot2)
library(data.table)
library(magrittr)
library(MetBrewer)

##

## Loading data
data_path = "~/Datascience/spatial_rnaseq/data/ND_Rexach_2024"
sc_object <- SingleCells(dir_data = data_path)
params_sc_min_quality()
sc_object <- bixverse::load_h5ad(
  object = sc_object,
  h5_path = file.path(data_path, "ND_Rexach_2024.h5ad")
)

var <- get_sc_var(sc_object)
setnames_sc(
  object = sc_object,
  table = "var",
  old = "feature_name",
  new = "gene_symbol"
)

var <- get_sc_var(sc_object)
var
ensembl_to_symbol <- setNames(var$gene_symbol, var$gene_id)
symbol_to_ensembl <- setNames(var$gene_id, var$gene_symbol)

## Quality control
## Gene set proportion

gs_of_interest <- list(
  MT = var[grepl("^MT-", gene_symbol), gene_id],
  Ribo = var[grepl("^RPS|^RPL", gene_symbol), gene_id]
)

sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  streaming = FALSE,
  .verbose = TRUE
)

sc_object[[1:5L]]

sc_object[["condition"]] <- dplyr::case_when(
  sc_object[[]]$disease == "Alzheimer disease" ~ "AD",
  sc_object[[]]$disease == "Pick disease" ~ "Pick",
  sc_object[[]]$disease == "progressive supranuclear palsy" ~ "PSP",
  sc_object[[]]$disease == "normal" ~ "CTRL"
)

## MAD outlier detection
qc_df <- sc_object[[c("cell_id", "donor_id", "lib_size", "nnz", "MT")]]

metrics <- list(
  log10_lib_size = log10(qc_df$lib_size),
  log10_nnz = log10(qc_df$nnz),
  MT = qc_df$MT
)

directions <- c(
  log10_lib_size = "twosided",
  log10_nnz = "twosided",
  MT = "above"
)

qc <- run_cell_qc(
  metrics = metrics,
  cells_to_keep = get_cells_to_keep(sc_object),
  directions = directions,
  threshold = 3,
  groups = qc_df$donor_id
)


## get_sc_obs --> the same naming?
qc_df <- get_data(qc)

plot_qc_violin(
  df = qc_df,
  grouping_column = "grp",
  variable = "log10_lib_size",
  log_scale = FALSE,
  outlier_column = "global_outlier",
  show_outlier = TRUE
)

plot_qc_density(
  df = qc_df,
  grouping_column = "grp",
  variable = "log10_lib_size",
  log_scale = FALSE
)

plot_joined_qc(
  df = qc_df,
  library_size = "log10_lib_size",
  nb_features = "log10_nnz",
  log_scale = FALSE
)
