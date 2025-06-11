library(devtools)
library(ggplot2)
library(magrittr)

rextendr::document()
# devtools::document()
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
  gtex_subgrp == "Brain - Caudate (basal ganglia)",
  sample_id
]
samples_to_keep_2 <- new_meta_data[
  gtex_subgrp == "Brain - Frontal Cortex (BA9)",
  sample_id
]
data_1 <- t(d)[samples_to_keep, ]
data_2 <- t(d)[samples_to_keep_2, ]
meta_data_1 <- new_meta_data[gtex_subgrp == "Brain - Caudate (basal ganglia)"]
meta_data_2 <- new_meta_data[gtex_subgrp == "Brain - Frontal Cortex (BA9)"]

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

devtools::load_all()

rextendr::document()

tictoc::tic()
cor_test <- bulk_coexp(raw_data = data_1, meta_data = meta_data_1) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  cor_module_processing(., cor_method = "spearman")

plot_hvgs(cor_test)

cor_test <- cor_module_check_epsilon(cor_test, rbf_func = "gaussian")

plot_epsilon_res(cor_test)

cor_test <- cor_module_coremo_clustering(
  object = cor_test,
  coremo_params = params_coremo(epsilon = 2)
)

params_coremo(epsilon = 2)

cluster_data <- cor_module_coremo_stability(object = cor_test)

devtools::load_all()

syn_data <- synthetic_signal_matrix()

data <- t(syn_data$mat)

meta_data <- data.table::data.table(
  sample_id = rownames(data),
  case_control = syn_data$group
)

cor_test <- bulk_coexp(raw_data = data, meta_data = meta_data) %>%
  preprocess_bulk_coexp(., hvg = 0.5) %>%
  cor_module_processing(., cor_method = "spearman") %>%
  cor_module_check_epsilon(., rbf_func = "gaussian")


dim(cor_test@processed_data$correlation_res$get_cor_matrix(.verbose = TRUE))

cor_test <- cor_module_coremo_clustering(
  object = cor_test,
  coremo_params = params_coremo(epsilon = 1)
)

saveRDS(
  cor_test@processed_data$feature_meta,
  "~/repos/shared/bixverse/inst/tinytest/test_data/hvg_data.rds"
)

cor_test <- cor_module_coremo_stability(
  object = cor_test
)
