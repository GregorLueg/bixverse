library(devtools)
library(ggplot2)
library(magrittr)

rextendr::document()
devtools::document()
devtools::load_all()

# Test on real data ------------------------------------------------------------

gtex_brain <- recount3::create_rse_manual(
  project = "BRAIN",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v29",
  type = "gene"
)

coldata <- SummarizedExperiment::colData(gtex_brain) |> as.data.frame()
rowdata <- SummarizedExperiment::rowData(gtex_brain) |> as.data.frame()

d <- edgeR::DGEList(SummarizedExperiment::assay(gtex_brain))
d <- edgeR::calcNormFactors(d, method = "upperquartile")
to_keep <- suppressWarnings(edgeR::filterByExpr(d))
d <- d[to_keep, ]
d <- edgeR::cpm(d, log = TRUE)

d <- as.matrix(d)

new_meta_data <- data.table::data.table(
  sample_id = rownames(coldata),
  case_control = "case",
  gtex_subgrp = coldata$gtex.smtsd
)

samples_to_keep <- new_meta_data[
  gtex_subgrp == "Brain - Hippocampus",
  sample_id
]
samples_to_keep_2 <- new_meta_data[gtex_subgrp == "Brain - Amygdala", sample_id]
data_1 <- t(d)[samples_to_keep, ]
data_2 <- t(d)[samples_to_keep_2, ]
meta_data_1 <- new_meta_data[gtex_subgrp == "Brain - Hippocampus"]
meta_data_2 <- new_meta_data[gtex_subgrp == "Brain - Amygdala"]

cor_test <- bulk_coexp(raw_data = data_1, meta_data = meta_data_1) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  cor_module_processing(., cor_method = "spearman")

devtools::document()

cor_test <- cor_module_check_epsilon(cor_test, rbf_func = "gaussian")

plot_epsilon_res(cor_test)

options(future.globals.maxSize = 2000 * 1024^2)
cor_test <- cor_module_check_res(
  cor_test,
  graph_params = params_cor_graph(epsilon = 1.5)
)

plot_resolution_res(cor_test)

devtools::load_all()

cor_test <- cor_module_final_modules(cor_test)

cor_test <- bulk_coexp(raw_data = data_2, meta_data = meta_data_2) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  cor_module_processing(., cor_method = "spearman")

# Differential correlation -----------------------------------------------------

devtools::load_all()
devtools::document()

dim(data_1)

tictoc::tic()
cor_test_2 <- bulk_coexp(raw_data = data_1, meta_data = meta_data_1) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  diffcor_module_processing(
    .,
    background_mat = data_2,
    cor_method = "spearman"
  ) %>%
  cor_module_check_res(.) %>%
  cor_module_final_modules(.)

cor_test_2 <- cor_module_final_modules(cor_test_2)
tictoc::toc()

# CoReMo -----------------------------------------------------------------------

cor_test <- bulk_coexp(raw_data = data_1, meta_data = meta_data_1) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  cor_module_processing(., cor_method = "spearman")

cor_test <- cor_module_check_epsilon(cor_test, rbf_func = "gaussian")

plot_epsilon_res(cor_test)

object = cor_test
epsilon = 2.5
rbf_func = "gaussian"
cluster_method = "ward.D"
cor_method = "spearman"
k_min = 1L
k_max = 100L
min_size = NULL
.verbose = FALSE

cor_res <- S7::prop(object, "processed_data")$correlation_res
cor_mat <- cor_res$get_cor_matrix()
dist_mat <- 1 - abs(cor_mat)
aff_mat <- rs_rbf_function_mat(
  x = dist_mat,
  epsilon = epsilon,
  rbf_type = rbf_func
)
dist_mat <- 1 - aff_mat

tree <- stats::hclust(as.dist(dist_mat), method = cluster_method)

# devtools::load_all()

tictoc::tic()
cutOpt <- tree_cut_iter(
  tree = tree,
  cor_mat = cor_mat,
  dist_mat = dist_mat,
  k_min = k_min,
  k_max = k_max,
  min_size = min_size,
  cor_method = cor_method
)
tictoc::toc()

k = 25L
tictoc::tic()
modules <-
  coremo_tree_cut(
    tree = tree,
    k = k,
    min_size = min_size,
    dist_mat = dist_mat,
    cor_method = cor_method
  ) %>%
  `names<-`(rownames(cor_mat))
tictoc::toc()

cluster_list <- split(names(modules), modules)

# rextendr::document()

tictoc::tic()
qc <- coremo_cluster_quality(modules = modules, cor_mat = cor_mat)
tictoc::toc()

tictoc::tic()
qc_2 <- rs_coremo_quality(
  cluster_genes = cluster_list,
  cor_mat = cor_mat,
  row_names = rownames(cor_mat),
  seed = 10101L
)
tictoc::toc()


coremo_cluster_quality_v2 <- function(modules, cor_mat, random_seed = 10101L) {
  # Checks
  checkmate::qassert(modules, c("S+", "I+"))
  checkmate::assertNamed(modules)
  checkmate::assertMatrix(cor_mat, mode = "numeric")
  checkmate::qassert(random_seed, "I1")
  # Function body
  cluster_list <- split(names(modules), modules)

  n_clusters <- length(cluster_list)
  result_list <- vector("list", n_clusters)

  for (i in seq_len(n_clusters)) {
    genes <- cluster_list[[i]]
    n <- length(genes)
    if (n < 2) {
      result_list[[i]] <- data.table::data.table(
        size = n,
        r2med = 1,
        r2mad = 0
      )
    } else {
      # Sample genes if needed
      scluster <- if (n > 1000) {
        sample(genes, 1000, replace = FALSE)
      } else {
        genes
      }

      # Extract submatrix more efficiently
      sub_cor_sq <- cor_mat[scluster, scluster, drop = FALSE]^2

      r2_values <- rs_coremo_quality(sub_cor_sq)

      result_list[[i]] <- data.table::data.table(
        size = n,
        r2med = r2_values$median,
        r2mad = r2_values$mad
      )
    }
  }

  data.table::rbindlist(result_list)
}
